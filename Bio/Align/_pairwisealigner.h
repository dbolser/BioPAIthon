/* Copyright 2018-2019 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */


#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "../_freethreading.h"

#define HORIZONTAL 0x1
#define VERTICAL 0x2
#define DIAGONAL 0x4

typedef enum {NeedlemanWunschSmithWaterman,
              Gotoh,
              WatermanSmithBeyer,
              FOGSAA,
              Unknown} Algorithm;

typedef enum {Global, Local, FOGSAA_Mode} Mode;

typedef struct {
    PyObject_HEAD
    Mode mode;
    Algorithm algorithm;
    double match;
    double mismatch;
    double epsilon;
    double open_internal_insertion_score;
    double extend_internal_insertion_score;
    double open_left_insertion_score;
    double extend_left_insertion_score;
    double open_right_insertion_score;
    double extend_right_insertion_score;
    double open_internal_deletion_score;
    double extend_internal_deletion_score;
    double open_left_deletion_score;
    double extend_left_deletion_score;
    double open_right_deletion_score;
    double extend_right_deletion_score;
    bool open_internal_insertion_score_set;
    bool extend_internal_insertion_score_set;
    bool open_left_insertion_score_set;
    bool extend_left_insertion_score_set;
    bool open_right_insertion_score_set;
    bool extend_right_insertion_score_set;
    bool open_internal_deletion_score_set;
    bool extend_internal_deletion_score_set;
    bool open_left_deletion_score_set;
    bool extend_left_deletion_score_set;
    bool open_right_deletion_score_set;
    bool extend_right_deletion_score_set;
    PyObject* insertion_score_function;
    PyObject* deletion_score_function;
    Py_buffer substitution_matrix;
    PyObject* alphabet;
    int wildcard;
} Aligner;


/* Release the references held by a snapshot taken by Aligner_snapshot.
 * Safe to call on a snapshot whose creation failed. */
static inline void
Aligner_snapshot_release(Aligner* snap)
{
    Py_CLEAR(snap->insertion_score_function);
    Py_CLEAR(snap->deletion_score_function);
    PyBuffer_Release(&snap->substitution_matrix);  /* no-op if obj is NULL */
}

/* Copy the aligner's settings into snap, for the duration of one call.
 *
 * A call can run Python code partway through (a gap function, or the
 * __repr__ of one), and that code may reconfigure the aligner: replace
 * the gap functions or the substitution matrix, freeing the old ones.
 * The snapshot therefore holds its own references to both gap functions
 * and its own buffer export of the substitution matrix, so the call
 * finishes with the settings it started with.
 *
 * The snapshot is a C struct, not a Python object: its PyObject_HEAD is a
 * stale copy of the aligner's.  Code given a snapshot as self (the
 * alignment kernels) must use it only as a C struct, and never pass it to
 * the Python API.  alphabet is set to NULL, as no reference to it is taken.
 *
 * Returns 0 on success.  On failure returns -1 with an exception set, and
 * snap holds no references, so Aligner_snapshot_release is still safe. */
static inline int
Aligner_snapshot(Aligner* self, Aligner* snap)
{
    PyObject* matrix;
    Py_buffer* view = &snap->substitution_matrix;
    int status;

    /* Copy the settings and take the references under the aligner's lock,
     * so that a setter in another thread cannot free an object between
     * the copy and the reference.  Nothing here calls into Python. */
    Py_BEGIN_CRITICAL_SECTION(self);
    *snap = *self;
    Py_XINCREF(snap->insertion_score_function);
    Py_XINCREF(snap->deletion_score_function);
    matrix = self->substitution_matrix.obj;
    Py_XINCREF(matrix);
    Py_END_CRITICAL_SECTION();
    snap->alphabet = NULL;
    view->obj = NULL;
    view->buf = NULL;
    if (!matrix) return 0;

    /* The exporter may run arbitrary code, so export the matrix outside
     * the lock, through the reference taken above. */
    status = PyObject_GetBuffer(matrix, view, PyBUF_FORMAT | PyBUF_ND);
    Py_DECREF(matrix);
    if (status < 0) {
        Aligner_snapshot_release(snap);
        return -1;
    }
    /* The setter checked the matrix when it was assigned, but its shape
     * or dtype may have been changed in place since then. */
    if (view->ndim != 2
     || view->shape[0] != view->shape[1]
     || view->len == 0
     || view->itemsize != sizeof(double)
     || view->format == NULL
     || strcmp(view->format, "d") != 0) {
        PyErr_SetString(PyExc_ValueError,
                        "substitution matrix is no longer a square matrix "
                        "of float values");
        Aligner_snapshot_release(snap);
        return -1;
    }
    return 0;
}
