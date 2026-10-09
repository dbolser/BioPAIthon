/* This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * Free-threading support shared by the C extensions.
 *
 * Every extension uses single-phase initialisation, so a module declares
 * that it can run without the GIL by calling Bio_module_gil_not_used() on
 * the module object in its PyInit_* function. This is CPython's documented
 * route for single-phase modules (PyUnstable_Module_SetGIL). Only call it
 * for a module that is safe without the GIL: a free-threaded interpreter
 * otherwise re-enables the GIL when the module is imported.
 */

#ifndef BIO_FREETHREADING_H
#define BIO_FREETHREADING_H

#include "Python.h"

/* Critical sections arrived in Python 3.13, and PyPy does not provide them.
 * Where they are missing there is always a GIL, so a critical section
 * reduces to the block it opens, which is exactly what CPython's own
 * definitions expand to on a GIL build. Expanding to braces keeps a
 * mismatched BEGIN/END a compile error on every build. */
#ifndef Py_BEGIN_CRITICAL_SECTION
#define Py_BEGIN_CRITICAL_SECTION(op) {
#define Py_END_CRITICAL_SECTION() }
#define Py_BEGIN_CRITICAL_SECTION2(a, b) {
#define Py_END_CRITICAL_SECTION2() }
#endif

/* Declare that module can run without the GIL. This does nothing unless
 * Py_GIL_DISABLED is defined, that is, unless this is a free-threaded build.
 * Returns 0 on success, or -1 with an exception set. */
static inline int
Bio_module_gil_not_used(PyObject *module)
{
#ifdef Py_GIL_DISABLED
    return PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#else
    (void)module;
    return 0;
#endif
}

#endif /* BIO_FREETHREADING_H */
