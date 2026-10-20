/* This code is part of the BioPAIthon distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */

/* Linear-space traceback for PairwiseAligner.align().
 *
 * Fork-owned, and #included once into _pairwisealigner.c after the
 * alignment kernels, so that upstream syncs only meet four hooks there.
 *
 * Above a threshold of traceback matrix size, align() returns a LinearPaths
 * object instead of building the matrix.  It yields the same paths in the
 * same order:
 *
 * - align() runs one forward pass with the kernel's per-cell arithmetic, so
 *   the score is bit-identical, saving every stride-th row as a checkpoint.
 * - The first path is recovered strip by strip, from the strip where it
 *   ends: recompute a strip's rows from its checkpoint, keeping their trace
 *   bits, and follow the traceback up through it.  Only the columns left
 *   of where the path leaves the strip are needed.  A strip over the block
 *   budget is bisected where that saves memory: sweep to its middle row,
 *   solve the bottom half, then the top half.  Same per-cell code, same
 *   trace bits, same path.  (Hirschberg's forward-and-reverse midpoints
 *   would pick another co-optimal path.)
 * - len() counts the paths in one more forward pass, recomputing each row's
 *   trace bits and keeping a row of counts.
 * - Later paths build the full matrix with the existing kernel and delegate
 *   to its PathGenerator, after checking it agrees.  So does the first
 *   path where the full matrix's own stops short of the origin, as Gotoh's
 *   can with infinite or NaN scores, or where its clean-up of a Gotoh local
 *   alignment removed a trace the path would take.
 *
 * It is used only where it holds less memory than the matrix.  A row of
 * trace bits costs less than a row of doubles, so a short target against
 * a long query has too few rows to save.
 *
 * Settings, substitution matrix and sequences are copied at align(), so
 * later changes to them cannot change the results.  Each algorithm
 * supplies a LinearKernel to the generic strip driver; Needleman-Wunsch,
 * Smith-Waterman and Gotoh global and local alignments have one.  The
 * threshold is 512 MiB by default; the private _set_traceback_limits
 * changes it, for testing.
 */


typedef struct LinearPaths LinearPaths;
typedef struct LinearWalker LinearWalker;

/* Return values of LinearKernel.walk and of linear_solve. */
#define LINEAR_EXITED 0   /* the walk reached the top row of the block */
#define LINEAR_ENDED 1    /* the path is complete */
#define LINEAR_STOPPED 2  /* the path stops before the origin, as the full
                           * matrix's does where a trace is empty */

typedef struct {
    /* A saved row holds ncarried scalars, then ndoubles per column. */
    size_t ncarried;
    size_t ndoubles;
    size_t trace_bytes;  /* bytes of trace bits per cell */
    size_t ncounts;      /* path counts per column, for len() */
    /* The full traceback matrix, as _matrix_size counts it: bytes per cell,
     * and row pointers per row. */
    size_t matrix_cell_bytes;
    size_t matrix_row_pointers;
    /* Fill row 0 for all columns, and its path counts if counts is not
     * NULL. */
    void (*start)(const LinearPaths* self, double* row, Py_ssize_t* counts);
    /* Compute rows r0+1 to r1 over columns 0 to c.  On entry row holds
     * row r0, on return row r1.  The trace bits of row r0+1+k are written
     * at T + k * tstride; tstride 0 overwrites a single scratch row.  If
     * counts is not NULL, c is the last column, and the path counts of
     * row r0 are updated to those of row r1 as in PathGenerator_length.
     * Releases the GIL.  Returns 0, or -1 if interrupted by a signal. */
    int (*sweep)(const LinearPaths* self, double* row, int r0, int r1, int c,
                 unsigned char* T, size_t tstride, Py_ssize_t* counts);
    /* The score of the alignment, from the last row.  Also sets end, 0 if
     * the full matrix has no paths, and the cell the first path ends in. */
    double (*score)(LinearPaths* self, const double* row);
    /* The number of paths, from the last row and its path counts, or
     * OVERFLOW_ERROR. */
    Py_ssize_t (*total)(const LinearPaths* self, const double* row,
                        const Py_ssize_t* counts);
    /* Follow the traceback up through the block holding the trace bits of
     * rows r0+1 to the walker's row, width cells per row.  Returns
     * LINEAR_EXITED with the walker in row r0, LINEAR_ENDED,
     * LINEAR_STOPPED, or -1 with an exception set. */
    int (*walk)(const LinearPaths* self, LinearWalker* walker,
                const unsigned char* T, int r0, size_t width);
    /* The kernel building the full traceback matrix. */
    PyObject* (*align)(Aligner* self, const int* sA, int nA,
                       const int* sB, int nB, unsigned char strand);
} LinearKernel;

/* What SMITHWATERMAN_ALIGN and GOTOH_LOCAL_ALIGN keep while they fill the
 * matrix, carried from row to row at the start of a saved row.  They mark
 * a cell as an end point if it is not a start point, its score is less
 * than epsilon below the highest score so far, and, for Smith-Waterman,
 * its trace has a diagonal move (every move into Gotoh's M is).
 * Where the score is more than epsilon above, a rise, they unmark every
 * cell from the previous rise on.  So the end points are the cells marked
 * before the first rise or from the last one on.  Cells are (i, j),
 * compared in row-major order; i = 0 is none. */
typedef struct {
    double maximum;      /* the highest score so far; for Gotoh, the score
                          * of the last rise */
    Py_ssize_t total;    /* for len(), the paths ending in the cells so far */
    int first[2];        /* the first marked cell */
    int rise[2];         /* the first rise */
    int last[2];         /* the last rise */
} LinearLocal;

#define LINEAR_LOCAL_DOUBLES \
    ((sizeof(LinearLocal) + sizeof(double) - 1) / sizeof(double))

/* How a LinearPaths object runs, fixed at align() from the limits. */
typedef struct {
    size_t rowsize;      /* doubles in a saved row */
    size_t stride;       /* rows between checkpoints */
    size_t count;        /* number of checkpoints */
    size_t depth;        /* most levels of bisection, each holding a row */
    size_t cells;        /* trace cells in a block */
    double peak;         /* most bytes held at once, an estimate */
} LinearPlan;

struct LinearPaths {
    PyObject_HEAD
    Aligner aligner;     /* snapshot of the settings, used as a C struct only */
    Py_ssize_t shape[2]; /* shape of the private copy of the substitution matrix */
    const LinearKernel* kernel;
    int* sA;
    int* sB;
    int nA;
    int nB;
    unsigned char strand;
    double score;
    int end;             /* the walker's trace at the end of the first path
                          * (see LinearWalker); 0 if there are no paths */
    int end_i;           /* the cell the first path ends in */
    int end_j;
    LinearLocal local;   /* local alignments: as the forward pass left it */
    LinearPlan plan;
    double* checkpoints; /* rows 0, stride, 2 * stride, ...; freed once the
                          * first path is known */
    PyObject* first;     /* the first path, once known */
    Py_ssize_t length;   /* the number of paths, once counted, or
                          * OVERFLOW_ERROR; 0 if not counted */
    PathGenerator* paths;/* the full traceback matrix, once built */
    bool primed;         /* paths has already yielded the first path */
    Py_ssize_t position; /* paths returned since the last reset */
    PyThread_type_lock lock;
    bool locked;         /* locked and owner are read and set in a critical
                          * section on the object */
    unsigned long owner;
};

struct LinearWalker {
    int i;               /* current cell */
    int j;
    int last;            /* direction of the move out of (i, j); 0 at the end */
    int trace;           /* for Gotoh, the matrices the path can be in at
                          * (i, j): the trace bits of the next cell on the
                          * path, or LinearPaths.end (M for a local
                          * alignment) at the end; for
                          * Smith-Waterman, the moves it may make out of
                          * (i, j): only DIAGONAL out of the end */
    int* points;         /* (i, j) where the path turns, last one first */
    size_t npoints;
    size_t capacity;
};

/* Limits, set by _set_traceback_limits; align() reads them under the GIL.
 * A negative threshold routes every alignment, for testing. */
static bool linear_enabled = true;
static Py_ssize_t linear_threshold = (Py_ssize_t)512 << 20;
static size_t linear_checkpoint_bytes = (size_t)32 << 20;
static size_t linear_block_bytes = (size_t)16 << 20;

#define LINEAR_INTERNAL_ERROR \
    PyErr_SetString(PyExc_RuntimeError, \
                    "internal error in the linear-space traceback; " \
                    "please report")

/* next(), len() and reset() hold the lock throughout, so that calls from
 * other threads wait while the GIL is released.  Signal handlers run while
 * it is held, so a call from the thread holding it fails rather than
 * waiting for itself.  Without the GIL, another thread could read locked
 * and owner while the thread taking the lock has set only one of them, and
 * so take itself for the holder; the critical sections make each pair of
 * reads and writes one step. */
static int
linear_lock(LinearPaths* self)
{
    const unsigned long thread = PyThread_get_thread_ident();
    bool held;
    Py_BEGIN_CRITICAL_SECTION(self);
    held = self->locked && self->owner == thread;
    Py_END_CRITICAL_SECTION();
    if (held) {
        PyErr_SetString(PyExc_RuntimeError,
                        "alignments used while they are being computed");
        return -1;
    }
    if (!PyThread_acquire_lock(self->lock, 0)) {
        Py_BEGIN_ALLOW_THREADS
        PyThread_acquire_lock(self->lock, 1);
        Py_END_ALLOW_THREADS
    }
    Py_BEGIN_CRITICAL_SECTION(self);
    self->locked = true;
    self->owner = thread;
    Py_END_CRITICAL_SECTION();
    return 0;
}

static void
linear_unlock(LinearPaths* self)
{
    Py_BEGIN_CRITICAL_SECTION(self);
    self->locked = false;
    Py_END_CRITICAL_SECTION();
    PyThread_release_lock(self->lock);
}


/* ---------------------- the walker ---------------------- */

/* Record the cell the walker is in. */
static int
linear_walker_add(LinearWalker* w)
{
    if (w->npoints == w->capacity) {
        const size_t capacity = w->capacity ? 2 * w->capacity : 64;
        int* points;
        if (capacity > PY_SSIZE_T_MAX / (2 * sizeof(int))) {
            PyErr_NoMemory();
            return -1;
        }
        points = PyMem_Realloc(w->points, capacity * 2 * sizeof(int));
        if (!points) {
            PyErr_NoMemory();
            return -1;
        }
        w->points = points;
        w->capacity = capacity;
    }
    w->points[2 * w->npoints] = w->i;
    w->points[2 * w->npoints + 1] = w->j;
    w->npoints++;
    return 0;
}

/* The cell the walker is in is entered by a move in direction d.  Record
 * the cell if the path turns there, as PathGenerator_create_path does. */
static int
linear_walker_step(LinearWalker* w, int d)
{
    if (d == w->last) return 0;
    w->last = d;
    return linear_walker_add(w);
}

/* The path as PathGenerator_create_path returns it: the coordinates at
 * every change of direction, in a tuple (target row, query row). */
static PyObject*
linear_walker_path(const LinearPaths* self, const LinearWalker* w)
{
    const Py_ssize_t n = (Py_ssize_t)w->npoints;
    const int nB = self->nB;
    Py_ssize_t k;
    PyObject* target_row;
    PyObject* query_row;
    PyObject* tuple = PyTuple_New(2);
    if (!tuple) return NULL;
    target_row = PyTuple_New(n);
    if (!target_row) goto error;
    PyTuple_SET_ITEM(tuple, 0, target_row);
    query_row = PyTuple_New(n);
    if (!query_row) goto error;
    PyTuple_SET_ITEM(tuple, 1, query_row);
    for (k = 0; k < n; k++) {
        const int* point = w->points + 2 * (n - 1 - k);
        PyObject* value = PyLong_FromLong(point[0]);
        if (!value) goto error;
        PyTuple_SET_ITEM(target_row, k, value);
        value = PyLong_FromLong(self->strand == '-' ? nB - point[1] : point[1]);
        if (!value) goto error;
        PyTuple_SET_ITEM(query_row, k, value);
    }
    return tuple;
error:
    Py_DECREF(tuple);
    return NULL;
}


/* ------------------- Needleman-Wunsch ------------------- */

/* SELECT_TRACE_NEEDLEMAN_WUNSCH, storing the trace bits in a row of
 * bytes, tr, rather than in M[i][j].  The arithmetic must stay
 * identical to it, and so must the order of the stores, or the trace
 * bits, and with them the paths, could differ from the full matrix. */
#define LINEAR_NW_CELL(hgap, vgap, align_score) \
    score = temp + (align_score); \
    trace = DIAGONAL; \
    temp = row[j-1] + hgap; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = HORIZONTAL; \
    } \
    else if (temp > score - epsilon) trace |= HORIZONTAL; \
    temp = row[j] + vgap; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = VERTICAL; \
    } \
    else if (temp > score - epsilon) trace |= VERTICAL; \
    temp = row[j]; \
    row[j] = score; \
    tr[j] = (unsigned char)trace;

/* One row of PathGenerator_needlemanwunsch_length, over the row's trace
 * bits in tr: the overflows, and the order of the additions, must stay as
 * there for len() to agree with it. */
static void
linear_nw_count(const unsigned char* tr, int nB, Py_ssize_t* counts)
{
    int j;
    int trace;
    Py_ssize_t term;
    Py_ssize_t count;
    Py_ssize_t temp;
    trace = tr[0];
    count = 0;
    if (trace & VERTICAL) SAFE_ADD(counts[0], count);
    temp = counts[0];
    counts[0] = count;
    for (j = 1; j <= nB; j++) {
        trace = tr[j];
        count = 0;
        if (trace & HORIZONTAL) SAFE_ADD(counts[j-1], count);
        if (trace & VERTICAL) SAFE_ADD(counts[j], count);
        if (trace & DIAGONAL) SAFE_ADD(temp, count);
        temp = counts[j];
        counts[j] = count;
    }
}

/* The rows of NEEDLEMANWUNSCH_ALIGN: the end gap scores swapped for the
 * reverse strand, the right horizontal gap score in the last row, the
 * right vertical one in the last column, and column 0 stored before the
 * row is computed, as there. */
#define LINEAR_NW_SWEEP(align_score) \
    const Aligner* const self = &lp->aligner; \
    const int nA = lp->nA; \
    const int nB = lp->nB; \
    const int* const sA = lp->sA; \
    const int* const sB = lp->sB; \
    const double gap_extend_A = self->extend_internal_insertion_score; \
    const double gap_extend_B = self->extend_internal_deletion_score; \
    const double epsilon = self->epsilon; \
    const bool forward = (lp->strand == '+'); \
    const double right_gap_extend_A = forward ? \
        self->extend_right_insertion_score : self->extend_left_insertion_score; \
    const double left_gap_extend_B = forward ? \
        self->extend_left_deletion_score : self->extend_right_deletion_score; \
    const double right_gap_extend_B = forward ? \
        self->extend_right_deletion_score : self->extend_left_deletion_score; \
    const int jmax = (c < nB) ? c : nB - 1; \
    int i; \
    int j; \
    int kA; \
    int kB; \
    int trace; \
    int interrupted = 0; \
    double hgap; \
    double score; \
    double temp; \
    unsigned char* tr; \
    PAIRWISE_NOGIL_BEGIN \
    for (i = r0 + 1; i <= r1; i++) { \
        PAIRWISE_NOGIL_CHECK(i) \
        tr = T + (size_t)(i - r0 - 1) * tstride; \
        hgap = (i == nA) ? right_gap_extend_A : gap_extend_A; \
        temp = row[0]; \
        row[0] = i * left_gap_extend_B; \
        tr[0] = VERTICAL; \
        kA = sA[i-1]; \
        for (j = 1; j <= jmax; j++) { \
            kB = sB[j-1]; \
            LINEAR_NW_CELL(hgap, gap_extend_B, align_score); \
        } \
        if (c == nB) { \
            kB = sB[j-1]; \
            LINEAR_NW_CELL(hgap, right_gap_extend_B, align_score); \
        } \
        if (counts) linear_nw_count(tr, nB, counts); \
    } \
    PAIRWISE_NOGIL_END \
    return interrupted ? -1 : 0;

static int
linear_nw_sweep_compare(const LinearPaths* lp, double* row, int r0, int r1,
                        int c, unsigned char* T, size_t tstride,
                        Py_ssize_t* counts)
{
    const double match = lp->aligner.match;
    const double mismatch = lp->aligner.mismatch;
    const int wildcard = lp->aligner.wildcard;
    LINEAR_NW_SWEEP(COMPARE_SCORE);
}

static int
linear_nw_sweep_matrix(const LinearPaths* lp, double* row, int r0, int r1,
                       int c, unsigned char* T, size_t tstride,
                       Py_ssize_t* counts)
{
    const Py_ssize_t n = lp->aligner.substitution_matrix.shape[0];
    LINEAR_NW_SWEEP(MATRIX_SCORE_NOGIL);
}

/* Row 0 of NEEDLEMANWUNSCH_ALIGN.  Its trace bits, as
 * PathGenerator_create_NWSW sets them, lead horizontally to the origin,
 * so each cell has one path. */
static void
linear_nw_start(const LinearPaths* lp, double* row, Py_ssize_t* counts)
{
    const double left_gap_extend_A = (lp->strand == '+') ?
        lp->aligner.extend_left_insertion_score :
        lp->aligner.extend_right_insertion_score;
    const int nB = lp->nB;
    int j;
    row[0] = 0;
    for (j = 1; j <= nB; j++) row[j] = j * left_gap_extend_A;
    if (counts) for (j = 0; j <= nB; j++) counts[j] = 1;
}

static double
linear_nw_score(LinearPaths* self, const double* row)
{
    self->end = M_MATRIX;  /* its one matrix; NEEDLEMANWUNSCH_ALIGN keeps
                            * the last trace */
    self->end_i = self->nA;
    self->end_j = self->nB;
    return row[self->nB];
}

static Py_ssize_t
linear_nw_total(const LinearPaths* self, const double* Py_UNUSED(row),
                const Py_ssize_t* counts)
{
    return counts[self->nB];
}

/* Follow the traceback as PathGenerator_next_needlemanwunsch does,
 * preferring horizontal, then vertical, then diagonal moves.  Row 0 is
 * not in any block: its trace is horizontal back to the origin. */
static int
linear_nw_walk(const LinearPaths* self, LinearWalker* w,
               const unsigned char* T, int r0, size_t width)
{
    while (w->i > r0) {
        const int trace = T[(size_t)(w->i - r0 - 1) * width + (size_t)w->j];
        if (trace & HORIZONTAL) {
            if (linear_walker_step(w, HORIZONTAL) < 0) return -1;
            w->j--;
        }
        else if (trace & VERTICAL) {
            if (linear_walker_step(w, VERTICAL) < 0) return -1;
            w->i--;
        }
        else if (trace & DIAGONAL) {
            if (linear_walker_step(w, DIAGONAL) < 0) return -1;
            w->i--;
            w->j--;
        }
        else {
            LINEAR_INTERNAL_ERROR;
            return -1;
        }
    }
    if (w->i > 0) return LINEAR_EXITED;
    if (w->j > 0) {
        if (linear_walker_step(w, HORIZONTAL) < 0) return -1;
        w->j = 0;
    }
    /* the origin always starts the path */
    if (linear_walker_add(w) < 0) return -1;
    return LINEAR_ENDED;
}

static const LinearKernel linear_nw_compare = {
    .ncarried = 0,
    .ndoubles = 1,
    .trace_bytes = 1,
    .ncounts = 1,
    .matrix_cell_bytes = sizeof(Trace),
    .matrix_row_pointers = 1,
    .start = linear_nw_start,
    .sweep = linear_nw_sweep_compare,
    .score = linear_nw_score,
    .total = linear_nw_total,
    .walk = linear_nw_walk,
    .align = Aligner_needlemanwunsch_align_compare,
};

static const LinearKernel linear_nw_matrix = {
    .ncarried = 0,
    .ndoubles = 1,
    .trace_bytes = 1,
    .ncounts = 1,
    .matrix_cell_bytes = sizeof(Trace),
    .matrix_row_pointers = 1,
    .start = linear_nw_start,
    .sweep = linear_nw_sweep_matrix,
    .score = linear_nw_score,
    .total = linear_nw_total,
    .walk = linear_nw_walk,
    .align = Aligner_needlemanwunsch_align_matrix,
};


/* ----------------------- Gotoh global ----------------------- */

/* A row holds the scores of M, Ix and Iy, in that order, for each column,
 * so that its first columns are a prefix.  A cell's trace bits take two
 * bytes, as in the full matrix: the trace of M, then those of Ix and Iy as
 * TraceGapsGotoh holds them, Ix in the low four bits. */
#define LINEAR_GOTOH_IX(tr, j) ((tr)[2*(size_t)(j)+1] & 0xf)
#define LINEAR_GOTOH_IY(tr, j) ((tr)[2*(size_t)(j)+1] >> 4)

/* SELECT_TRACE_GOTOH_GLOBAL_ALIGN and SELECT_TRACE_GOTOH_GLOBAL_GAP,
 * without the store of the trace.  The arithmetic must stay identical to
 * them, or the trace bits, and with them the paths, could differ from the
 * full matrix. */
#define LINEAR_GOTOH_ALIGN \
    trace = M_MATRIX; \
    score = M_temp; \
    temp = Ix_temp; \
    if (temp > score + epsilon) { \
        score = Ix_temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = Iy_temp; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX;

#define LINEAR_GOTOH_GAP(score1, score2, score3) \
    trace = M_MATRIX; \
    score = score1; \
    temp = score2; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = score3; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX;

/* A cell of GOTOH_GLOBAL_ALIGN, with the gap scores of its row and column,
 * in the order of its statements.  M, Ix and Iy of column j are row[3*j],
 * row[3*j+1] and row[3*j+2].  Indices are in size_t, as 3*j can be over
 * INT_MAX. */
#define LINEAR_GOTOH_CELL(open_B, extend_B, open_A, extend_A, align_score) \
    LINEAR_GOTOH_ALIGN; \
    tr[2*(size_t)j] = (unsigned char)trace; \
    M_temp = row[3*(size_t)j]; \
    row[3*(size_t)j] = score + (align_score); \
    LINEAR_GOTOH_GAP(M_temp + open_B, \
                     row[3*(size_t)j+1] + extend_B, \
                     row[3*(size_t)j+2] + open_B); \
    tr[2*(size_t)j+1] = (unsigned char)trace; \
    Ix_temp = row[3*(size_t)j+1]; \
    row[3*(size_t)j+1] = score; \
    LINEAR_GOTOH_GAP(row[3*(size_t)j-3] + open_A, \
                     row[3*(size_t)j-2] + open_A, \
                     row[3*(size_t)j-1] + extend_A); \
    tr[2*(size_t)j+1] |= (unsigned char)(trace << 4); \
    Iy_temp = row[3*(size_t)j+2]; \
    row[3*(size_t)j+2] = score;

/* One row of PathGenerator_gotoh_global_length, over the row's trace bits
 * in tr, with the counts of M, Ix and Iy one after the other in counts:
 * the overflows, and the order of the additions, must stay as there for
 * len() to agree with it. */
static void
linear_gotoh_count(const unsigned char* tr, int nB, Py_ssize_t* counts)
{
    Py_ssize_t* const M_counts = counts;
    Py_ssize_t* const Ix_counts = M_counts + nB + 1;
    Py_ssize_t* const Iy_counts = Ix_counts + nB + 1;
    int j;
    int trace;
    Py_ssize_t count;
    Py_ssize_t term;
    Py_ssize_t M_temp;
    Py_ssize_t Ix_temp;
    Py_ssize_t Iy_temp;
    M_temp = M_counts[0];
    M_counts[0] = 0;
    Ix_temp = Ix_counts[0];
    Ix_counts[0] = 1;
    Iy_temp = Iy_counts[0];
    Iy_counts[0] = 0;
    for (j = 1; j <= nB; j++) {
        count = 0;
        trace = tr[2*(size_t)j];
        if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
        if (trace & Ix_MATRIX) SAFE_ADD(Ix_temp, count);
        if (trace & Iy_MATRIX) SAFE_ADD(Iy_temp, count);
        M_temp = M_counts[j];
        M_counts[j] = count;
        count = 0;
        trace = LINEAR_GOTOH_IX(tr, j);
        if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
        if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j], count);
        if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j], count);
        Ix_temp = Ix_counts[j];
        Ix_counts[j] = count;
        count = 0;
        trace = LINEAR_GOTOH_IY(tr, j);
        if (trace & M_MATRIX) SAFE_ADD(M_counts[j-1], count);
        if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j-1], count);
        if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j-1], count);
        Iy_temp = Iy_counts[j];
        Iy_counts[j] = count;
    }
}

/* The rows of GOTOH_GLOBAL_ALIGN: the end gap scores swapped for the
 * reverse strand, the right insertion scores in the last row, the right
 * deletion scores in the last column, and column 0 as there, with the
 * trace bits PathGenerator_create_Gotoh gives it. */
#define LINEAR_GOTOH_SWEEP(align_score) \
    const Aligner* const self = &lp->aligner; \
    const int nA = lp->nA; \
    const int nB = lp->nB; \
    const int* const sA = lp->sA; \
    const int* const sB = lp->sB; \
    const double gap_open_A = self->open_internal_insertion_score; \
    const double gap_open_B = self->open_internal_deletion_score; \
    const double gap_extend_A = self->extend_internal_insertion_score; \
    const double gap_extend_B = self->extend_internal_deletion_score; \
    const double epsilon = self->epsilon; \
    const bool forward = (lp->strand == '+'); \
    const double left_gap_open_B = forward ? \
        self->open_left_deletion_score : self->open_right_deletion_score; \
    const double left_gap_extend_B = forward ? \
        self->extend_left_deletion_score : self->extend_right_deletion_score; \
    const double right_gap_open_A = forward ? \
        self->open_right_insertion_score : self->open_left_insertion_score; \
    const double right_gap_extend_A = forward ? \
        self->extend_right_insertion_score : self->extend_left_insertion_score; \
    const double right_gap_open_B = forward ? \
        self->open_right_deletion_score : self->open_left_deletion_score; \
    const double right_gap_extend_B = forward ? \
        self->extend_right_deletion_score : self->extend_left_deletion_score; \
    const int jmax = (c < nB) ? c : nB - 1; \
    int i; \
    int j; \
    int kA; \
    int kB; \
    int trace; \
    int interrupted = 0; \
    double open_A; \
    double extend_A; \
    double score; \
    double temp; \
    double M_temp; \
    double Ix_temp; \
    double Iy_temp; \
    unsigned char* tr; \
    PAIRWISE_NOGIL_BEGIN \
    for (i = r0 + 1; i <= r1; i++) { \
        PAIRWISE_NOGIL_CHECK(i) \
        tr = T + (size_t)(i - r0 - 1) * tstride; \
        if (i == nA) { \
            open_A = right_gap_open_A; \
            extend_A = right_gap_extend_A; \
        } \
        else { \
            open_A = gap_open_A; \
            extend_A = gap_extend_A; \
        } \
        kA = sA[i-1]; \
        M_temp = row[0]; \
        Ix_temp = row[1]; \
        Iy_temp = row[2]; \
        row[0] = -DBL_MAX; \
        row[1] = left_gap_open_B + left_gap_extend_B * (i-1); \
        row[2] = -DBL_MAX; \
        tr[0] = 0; \
        tr[1] = (i == 1) ? M_MATRIX : Ix_MATRIX; \
        for (j = 1; j <= jmax; j++) { \
            kB = sB[j-1]; \
            LINEAR_GOTOH_CELL(gap_open_B, gap_extend_B, open_A, extend_A, \
                              align_score); \
        } \
        if (c == nB) { \
            kB = sB[j-1]; \
            LINEAR_GOTOH_CELL(right_gap_open_B, right_gap_extend_B, \
                              open_A, extend_A, align_score); \
        } \
        if (counts) linear_gotoh_count(tr, nB, counts); \
    } \
    PAIRWISE_NOGIL_END \
    return interrupted ? -1 : 0;

static int
linear_gotoh_sweep_compare(const LinearPaths* lp, double* row, int r0, int r1,
                           int c, unsigned char* T, size_t tstride,
                           Py_ssize_t* counts)
{
    const double match = lp->aligner.match;
    const double mismatch = lp->aligner.mismatch;
    const int wildcard = lp->aligner.wildcard;
    LINEAR_GOTOH_SWEEP(COMPARE_SCORE);
}

static int
linear_gotoh_sweep_matrix(const LinearPaths* lp, double* row, int r0, int r1,
                          int c, unsigned char* T, size_t tstride,
                          Py_ssize_t* counts)
{
    const Py_ssize_t n = lp->aligner.substitution_matrix.shape[0];
    LINEAR_GOTOH_SWEEP(MATRIX_SCORE_NOGIL);
}

/* Row 0 of GOTOH_GLOBAL_ALIGN, and its path counts as
 * PathGenerator_gotoh_global_length sets them. */
static void
linear_gotoh_start(const LinearPaths* lp, double* row, Py_ssize_t* counts)
{
    const bool forward = (lp->strand == '+');
    const double left_gap_open_A = forward ?
        lp->aligner.open_left_insertion_score :
        lp->aligner.open_right_insertion_score;
    const double left_gap_extend_A = forward ?
        lp->aligner.extend_left_insertion_score :
        lp->aligner.extend_right_insertion_score;
    const int nB = lp->nB;
    int j;
    row[0] = 0;
    row[1] = -DBL_MAX;
    row[2] = -DBL_MAX;
    for (j = 1; j <= nB; j++) {
        row[3*(size_t)j] = -DBL_MAX;
        row[3*(size_t)j+1] = -DBL_MAX;
        row[3*(size_t)j+2] = left_gap_open_A + left_gap_extend_A * (j-1);
    }
    if (counts) {
        Py_ssize_t* const M_counts = counts;
        Py_ssize_t* const Ix_counts = M_counts + nB + 1;
        Py_ssize_t* const Iy_counts = Ix_counts + nB + 1;
        M_counts[0] = 1;
        Ix_counts[0] = 0;
        Iy_counts[0] = 0;
        for (j = 1; j <= nB; j++) {
            M_counts[j] = 0;
            Ix_counts[j] = 0;
            Iy_counts[j] = 1;
        }
    }
}

/* The score as GOTOH_GLOBAL_ALIGN takes it from the last cell, and the
 * matrices whose trace there it keeps: those not more than epsilon below
 * the score.  Every trace of that cell has a bit set until then, so these
 * are the matrices PathGenerator_next_gotoh_global and
 * PathGenerator_gotoh_global_length find the paths ending in. */
static double
linear_gotoh_score(LinearPaths* self, const double* row)
{
    const double epsilon = self->aligner.epsilon;
    const double* const cell = row + 3 * (size_t)self->nB;
    double score;
    double temp;
    SELECT_SCORE_GLOBAL(cell[0], cell[1], cell[2]);
    self->end = 0;
    if (!(cell[0] < score - epsilon)) self->end |= M_MATRIX;
    if (!(cell[1] < score - epsilon)) self->end |= Ix_MATRIX;
    if (!(cell[2] < score - epsilon)) self->end |= Iy_MATRIX;
    self->end_i = self->nA;
    self->end_j = self->nB;
    return score;
}

static Py_ssize_t
linear_gotoh_total(const LinearPaths* self, const double* Py_UNUSED(row),
                   const Py_ssize_t* counts)
{
    const int nB = self->nB;
    const Py_ssize_t* const M_counts = counts;
    const Py_ssize_t* const Ix_counts = M_counts + nB + 1;
    const Py_ssize_t* const Iy_counts = Ix_counts + nB + 1;
    Py_ssize_t count = 0;
    Py_ssize_t term;
    if (self->end & M_MATRIX) SAFE_ADD(M_counts[nB], count);
    if (self->end & Ix_MATRIX) SAFE_ADD(Ix_counts[nB], count);
    if (self->end & Iy_MATRIX) SAFE_ADD(Iy_counts[nB], count);
    return count;
}

/* Follow the traceback as PathGenerator_next_gotoh_global does: in each
 * cell, take M if the walker's trace allows it, else Ix, else Iy, and read
 * that matrix's trace of the cell.  An empty trace ends the path.  That is
 * at the origin, or in column 0 or row 0 where the full matrix's path
 * stops short of the origin.  Row 0 is not in any block: its Iy traces
 * lead horizontally to the origin, its M and Ix traces are empty. */
static int
linear_gotoh_walk(const LinearPaths* Py_UNUSED(self), LinearWalker* w,
                  const unsigned char* T, int r0, size_t width)
{
    while (w->i > r0) {
        const unsigned char* tr = T + (size_t)(w->i - r0 - 1) * width * 2;
        const size_t j = (size_t)w->j;
        if (w->trace & M_MATRIX) {
            if (linear_walker_step(w, DIAGONAL) < 0) return -1;
            w->trace = tr[2*j];
            w->i--;
            w->j--;
        }
        else if (w->trace & Ix_MATRIX) {
            if (linear_walker_step(w, VERTICAL) < 0) return -1;
            w->trace = LINEAR_GOTOH_IX(tr, j);
            w->i--;
        }
        else if (w->trace & Iy_MATRIX) {
            if (linear_walker_step(w, HORIZONTAL) < 0) return -1;
            w->trace = LINEAR_GOTOH_IY(tr, j);
            w->j--;
        }
        else {
            LINEAR_INTERNAL_ERROR;
            return -1;
        }
        /* only M and Iy in column 0 have an empty trace in a block */
        if (!w->trace) return LINEAR_STOPPED;
    }
    if (w->i > 0) return LINEAR_EXITED;
    if (w->j > 0) {
        if (w->trace & (M_MATRIX | Ix_MATRIX)) return LINEAR_STOPPED;
        if (linear_walker_step(w, HORIZONTAL) < 0) return -1;
        w->j = 0;
    }
    /* the origin always starts the path */
    if (linear_walker_add(w) < 0) return -1;
    return LINEAR_ENDED;
}

static const LinearKernel linear_gotoh_compare = {
    .ncarried = 0,
    .ndoubles = 3,
    .trace_bytes = 2,
    .ncounts = 3,
    .matrix_cell_bytes = sizeof(Trace) + sizeof(TraceGapsGotoh),
    .matrix_row_pointers = 2,
    .start = linear_gotoh_start,
    .sweep = linear_gotoh_sweep_compare,
    .score = linear_gotoh_score,
    .total = linear_gotoh_total,
    .walk = linear_gotoh_walk,
    .align = Aligner_gotoh_global_align_compare,
};

static const LinearKernel linear_gotoh_matrix = {
    .ncarried = 0,
    .ndoubles = 3,
    .trace_bytes = 2,
    .ncounts = 3,
    .matrix_cell_bytes = sizeof(Trace) + sizeof(TraceGapsGotoh),
    .matrix_row_pointers = 2,
    .start = linear_gotoh_start,
    .sweep = linear_gotoh_sweep_matrix,
    .score = linear_gotoh_score,
    .total = linear_gotoh_total,
    .walk = linear_gotoh_walk,
    .align = Aligner_gotoh_global_align_matrix,
};

/* ----------------------- Smith-Waterman ----------------------- */

/* Whether cell (i, j) comes before cell in row-major order. */
static inline bool
linear_before(int i, int j, const int* cell)
{
    return i < cell[0] || (i == cell[0] && j < cell[1]);
}

/* The end point test of SELECT_TRACE_SMITH_WATERMAN_HVD and
 * SELECT_TRACE_SMITH_WATERMAN_D.  Where they unmark the cells since the
 * previous rise, this records the rise instead (see LinearLocal). */
#define LINEAR_SW_MARK \
    else if (trace & DIAGONAL && score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            if (!local.rise[0]) { \
                local.rise[0] = i; \
                local.rise[1] = j; \
            } \
            local.last[0] = i; \
            local.last[1] = j; \
        } \
        if (!local.first[0]) { \
            local.first[0] = i; \
            local.first[1] = j; \
        } \
        trace |= ENDPOINT; \
    }

/* SELECT_TRACE_SMITH_WATERMAN_HVD, storing the trace bits in a row of
 * bytes, tr, rather than in M[i][j].  The arithmetic must stay identical
 * to it, and so must the order of the stores, or the trace bits, and with
 * them the paths, could differ from the full matrix. */
#define LINEAR_SW_HVD(align_score) \
    trace = DIAGONAL; \
    score = temp + (align_score); \
    temp = row[j-1] + gap_extend_A; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = HORIZONTAL; \
    } \
    else if (temp > score - epsilon) trace |= HORIZONTAL; \
    temp = row[j] + gap_extend_B; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = VERTICAL; \
    } \
    else if (temp > score - epsilon) trace |= VERTICAL; \
    if (score < epsilon) { \
        score = 0; \
        trace = STARTPOINT; \
    } \
    LINEAR_SW_MARK \
    tr[j] = (unsigned char)trace; \
    if (score > maximum) maximum = score; \
    temp = row[j]; \
    row[j] = score;

/* SELECT_TRACE_SMITH_WATERMAN_D, likewise. */
#define LINEAR_SW_D(align_score) \
    score = temp + (align_score); \
    trace = DIAGONAL; \
    if (score < epsilon) { \
        score = 0; \
    } \
    LINEAR_SW_MARK \
    tr[j] = (unsigned char)trace; \
    if (score > maximum) maximum = score; \
    temp = row[j]; \
    row[j] = score;

/* One row of PathGenerator_smithwaterman_length, over the row's trace bits
 * in tr: the overflows, and the order of the additions, must stay as there
 * for len() to agree with it.  That counts over the traces as
 * SMITHWATERMAN_ALIGN leaves them, after two changes.  It unmarks the cells
 * marked from the first rise to before the last (see LinearLocal).  And so
 * that no path goes on past an end point, it removes each trace leading
 * back to one, and clears the cells then left without a trace, whose count
 * is 0.  Here an end point keeps a count of 0 instead, which adds the same
 * to every sum. */
static void
linear_sw_count(const LinearPaths* self, const unsigned char* tr, int i,
                Py_ssize_t* counts, Py_ssize_t* total)
{
    const int nB = self->nB;
    const int* const rise = self->local.rise;
    const int* const last = self->local.last;
    int j;
    int trace;
    Py_ssize_t term;
    Py_ssize_t count;
    Py_ssize_t temp;
    Py_ssize_t sum = *total;
    temp = counts[0];
    counts[0] = 1;
    for (j = 1; j <= nB; j++) {
        trace = tr[j];
        if (trace & ENDPOINT && rise[0] && !linear_before(i, j, rise)
                             && linear_before(i, j, last))
            trace &= ~ENDPOINT;
        count = 0;
        if (trace & DIAGONAL) SAFE_ADD(temp, count);
        if (trace & ENDPOINT) SAFE_ADD(count, sum);
        if (trace & HORIZONTAL) SAFE_ADD(counts[j-1], count);
        if (trace & VERTICAL) SAFE_ADD(counts[j], count);
        temp = counts[j];
        if (count == 0 && (trace & STARTPOINT)) count = 1;
        counts[j] = (trace & ENDPOINT) ? 0 : count;
    }
    *total = sum;
}

/* The rows of SMITHWATERMAN_ALIGN: only diagonal moves into the last row
 * and the last column, and column 0 start points, as
 * PathGenerator_create_NWSW makes them.  A saved row holds a LinearLocal,
 * then the scores. */
#define LINEAR_SW_SWEEP(align_score) \
    const Aligner* const self = &lp->aligner; \
    const int nA = lp->nA; \
    const int nB = lp->nB; \
    const int* const sA = lp->sA; \
    const int* const sB = lp->sB; \
    const double gap_extend_A = self->extend_internal_insertion_score; \
    const double gap_extend_B = self->extend_internal_deletion_score; \
    const double epsilon = self->epsilon; \
    double* const row = saved + LINEAR_LOCAL_DOUBLES; \
    LinearLocal local; \
    double maximum; \
    int i; \
    int j; \
    int jmax; \
    int kA; \
    int kB; \
    int trace; \
    int interrupted = 0; \
    double score; \
    double temp; \
    unsigned char* tr; \
    memcpy(&local, saved, sizeof(LinearLocal)); \
    maximum = local.maximum; \
    PAIRWISE_NOGIL_BEGIN \
    for (i = r0 + 1; i <= r1; i++) { \
        PAIRWISE_NOGIL_CHECK(i) \
        tr = T + (size_t)(i - r0 - 1) * tstride; \
        tr[0] = STARTPOINT; \
        kA = sA[i-1]; \
        temp = 0; \
        jmax = (i == nA) ? 0 : (c < nB) ? c : nB - 1; \
        for (j = 1; j <= jmax; j++) { \
            kB = sB[j-1]; \
            LINEAR_SW_HVD(align_score); \
        } \
        for ( ; j <= c; j++) { \
            kB = sB[j-1]; \
            LINEAR_SW_D(align_score); \
        } \
        if (counts) linear_sw_count(lp, tr, i, counts, &local.total); \
    } \
    PAIRWISE_NOGIL_END \
    local.maximum = maximum; \
    memcpy(saved, &local, sizeof(LinearLocal)); \
    return interrupted ? -1 : 0;

static int
linear_sw_sweep_compare(const LinearPaths* lp, double* saved, int r0, int r1,
                        int c, unsigned char* T, size_t tstride,
                        Py_ssize_t* counts)
{
    const double match = lp->aligner.match;
    const double mismatch = lp->aligner.mismatch;
    const int wildcard = lp->aligner.wildcard;
    LINEAR_SW_SWEEP(COMPARE_SCORE);
}

static int
linear_sw_sweep_matrix(const LinearPaths* lp, double* saved, int r0, int r1,
                       int c, unsigned char* T, size_t tstride,
                       Py_ssize_t* counts)
{
    const Py_ssize_t n = lp->aligner.substitution_matrix.shape[0];
    LINEAR_SW_SWEEP(MATRIX_SCORE_NOGIL);
}

/* Row 0 of SMITHWATERMAN_ALIGN, nothing kept yet, and the path counts of
 * row 0 as PathGenerator_smithwaterman_length sets them. */
static void
linear_sw_start(const LinearPaths* lp, double* saved, Py_ssize_t* counts)
{
    const int nB = lp->nB;
    double* const row = saved + LINEAR_LOCAL_DOUBLES;
    const LinearLocal local = {0};
    int j;
    memcpy(saved, &local, sizeof(LinearLocal));
    for (j = 0; j <= nB; j++) row[j] = 0;
    if (counts) for (j = 0; j <= nB; j++) counts[j] = 1;
}

/* The score as SMITHWATERMAN_ALIGN and GOTOH_LOCAL_ALIGN return it, the
 * highest, and the end of the first path, with the walker's trace end
 * there.  PathGenerator_next_smithwaterman and
 * PathGenerator_next_gotoh_local end it at the first end point in
 * row-major order: the first marked cell if it comes before the first
 * rise, else the last rise (see LinearLocal).  With no marked cell there
 * are no paths. */
static double
linear_local_score(LinearPaths* self, const double* saved, int end)
{
    LinearLocal* const local = &self->local;
    memcpy(local, saved, sizeof(LinearLocal));
    self->end = 0;
    if (local->first[0]) {
        const bool unmarked = local->first[0] == local->rise[0]
                           && local->first[1] == local->rise[1];
        const int* const cell = unmarked ? local->last : local->first;
        self->end = end;
        self->end_i = cell[0];
        self->end_j = cell[1];
    }
    return local->maximum;
}

/* Every cell has a trace, to cells before it in row-major order, and no
 * cell before the first end point is one.  So SMITHWATERMAN_ALIGN removes
 * no trace from the cells before it (see linear_sw_count), and the walk
 * can follow them as the sweep computes them.  A marked cell scores above
 * 0, so with a highest score of 0, where the full matrix's generator has
 * no paths, there are none. */
static double
linear_sw_score(LinearPaths* self, const double* saved)
{
    return linear_local_score(self, saved, DIAGONAL);
}

static Py_ssize_t
linear_local_total(const LinearPaths* Py_UNUSED(self), const double* saved,
                   const Py_ssize_t* Py_UNUSED(counts))
{
    LinearLocal local;
    memcpy(&local, saved, sizeof(LinearLocal));
    return local.total;
}

/* Follow the traceback as PathGenerator_next_smithwaterman does: out of the
 * end diagonally, as it excludes paths ending in a gap, then preferring
 * horizontal, then vertical, then diagonal moves, to a start point.  Row 0
 * is not in any block: its cells are start points. */
static int
linear_sw_walk(const LinearPaths* Py_UNUSED(self), LinearWalker* w,
               const unsigned char* T, int r0, size_t width)
{
    while (w->i > r0) {
        const int trace = w->trace
                        & T[(size_t)(w->i - r0 - 1) * width + (size_t)w->j];
        w->trace = HORIZONTAL | VERTICAL | DIAGONAL | STARTPOINT;
        if (trace & HORIZONTAL) {
            if (linear_walker_step(w, HORIZONTAL) < 0) return -1;
            w->j--;
        }
        else if (trace & VERTICAL) {
            if (linear_walker_step(w, VERTICAL) < 0) return -1;
            w->i--;
        }
        else if (trace & DIAGONAL) {
            if (linear_walker_step(w, DIAGONAL) < 0) return -1;
            w->i--;
            w->j--;
        }
        else if (trace & STARTPOINT) {
            if (linear_walker_add(w) < 0) return -1;
            return LINEAR_ENDED;
        }
        else {
            LINEAR_INTERNAL_ERROR;
            return -1;
        }
    }
    if (w->i > 0) return LINEAR_EXITED;
    if (linear_walker_add(w) < 0) return -1;
    return LINEAR_ENDED;
}

static const LinearKernel linear_sw_compare = {
    .ncarried = LINEAR_LOCAL_DOUBLES,
    .ndoubles = 1,
    .trace_bytes = 1,
    .ncounts = 1,
    .matrix_cell_bytes = sizeof(Trace),
    .matrix_row_pointers = 1,
    .start = linear_sw_start,
    .sweep = linear_sw_sweep_compare,
    .score = linear_sw_score,
    .total = linear_local_total,
    .walk = linear_sw_walk,
    .align = Aligner_smithwaterman_align_compare,
};

static const LinearKernel linear_sw_matrix = {
    .ncarried = LINEAR_LOCAL_DOUBLES,
    .ndoubles = 1,
    .trace_bytes = 1,
    .ncounts = 1,
    .matrix_cell_bytes = sizeof(Trace),
    .matrix_row_pointers = 1,
    .start = linear_sw_start,
    .sweep = linear_sw_sweep_matrix,
    .score = linear_sw_score,
    .total = linear_local_total,
    .walk = linear_sw_walk,
    .align = Aligner_smithwaterman_align_matrix,
};

/* ----------------------- Gotoh local ----------------------- */

/* SELECT_TRACE_GOTOH_LOCAL_ALIGN, storing the trace in tr rather than in
 * M[i][j].  Where it unmarks the cells since the previous rise, this
 * records the rise instead (see LinearLocal).  The arithmetic must stay
 * identical to it, or the trace bits, and with them the paths, could
 * differ from the full matrix. */
#define LINEAR_GOTOH_LOCAL_ALIGN(align_score) \
    trace = M_MATRIX; \
    score = M_temp; \
    if (Ix_temp > score + epsilon) { \
        score = Ix_temp; \
        trace = Ix_MATRIX; \
    } \
    else if (Ix_temp > score - epsilon) trace |= Ix_MATRIX; \
    if (Iy_temp > score + epsilon) { \
        score = Iy_temp; \
        trace = Iy_MATRIX; \
    } \
    else if (Iy_temp > score - epsilon) trace |= Iy_MATRIX; \
    score += (align_score); \
    if (score < epsilon) { \
        score = 0; \
        trace = STARTPOINT; \
    } \
    else if (score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            maximum = score; \
            if (!local.rise[0]) { \
                local.rise[0] = i; \
                local.rise[1] = j; \
            } \
            local.last[0] = i; \
            local.last[1] = j; \
        } \
        if (!local.first[0]) { \
            local.first[0] = i; \
            local.first[1] = j; \
        } \
        trace |= ENDPOINT; \
    } \
    tr[2*(size_t)j] = (unsigned char)trace;

/* SELECT_TRACE_GOTOH_LOCAL_GAP, without the store of the trace. */
#define LINEAR_GOTOH_LOCAL_GAP(score1, score2, score3) \
    LINEAR_GOTOH_GAP(score1, score2, score3) \
    if (score < epsilon) { \
        score = -DBL_MAX; \
        trace = 0; \
    }

/* A cell of GOTOH_LOCAL_ALIGN before the last row and column, in the order
 * of its statements, with M, Ix and Iy of column j at row[3*j], row[3*j+1]
 * and row[3*j+2]. */
#define LINEAR_GOTOH_LOCAL_CELL(align_score) \
    LINEAR_GOTOH_LOCAL_ALIGN(align_score) \
    M_temp = row[3*(size_t)j]; \
    row[3*(size_t)j] = score; \
    LINEAR_GOTOH_LOCAL_GAP(M_temp + gap_open_B, \
                           row[3*(size_t)j+1] + gap_extend_B, \
                           row[3*(size_t)j+2] + gap_open_B); \
    tr[2*(size_t)j+1] = (unsigned char)trace; \
    Ix_temp = row[3*(size_t)j+1]; \
    row[3*(size_t)j+1] = score; \
    LINEAR_GOTOH_LOCAL_GAP(row[3*(size_t)j-3] + gap_open_A, \
                           row[3*(size_t)j-2] + gap_open_A, \
                           row[3*(size_t)j-1] + gap_extend_A); \
    tr[2*(size_t)j+1] |= (unsigned char)(trace << 4); \
    Iy_temp = row[3*(size_t)j+2]; \
    row[3*(size_t)j+2] = score;

/* A cell of GOTOH_LOCAL_ALIGN in the last row or column, where Ix and Iy
 * are 0 without a trace. */
#define LINEAR_GOTOH_LOCAL_EDGE(align_score) \
    LINEAR_GOTOH_LOCAL_ALIGN(align_score) \
    M_temp = row[3*(size_t)j]; \
    row[3*(size_t)j] = score; \
    Ix_temp = row[3*(size_t)j+1]; \
    row[3*(size_t)j+1] = 0; \
    Iy_temp = row[3*(size_t)j+2]; \
    row[3*(size_t)j+2] = 0; \
    tr[2*(size_t)j+1] = 0;

/* One row of PathGenerator_gotoh_local_length, over the row's trace bits
 * in tr, with the counts of M, Ix and Iy one after the other in counts:
 * the overflows, and the order of the additions, must stay as there for
 * len() to agree with it.  That counts over the traces as GOTOH_LOCAL_ALIGN
 * leaves them, after two changes, as for Smith-Waterman (see
 * linear_sw_count).  It unmarks the cells marked from the first rise to
 * before the last.  And it removes each trace to an end point's M, or to a
 * matrix that then no start point reaches, which it leaves without a
 * trace.  Here an end point's M passes on a count of 0 instead; so such a
 * matrix counts 0 too, and the same is added to every sum. */
static void
linear_gotoh_local_count(const LinearPaths* self, const unsigned char* tr,
                         int i, Py_ssize_t* counts, Py_ssize_t* total)
{
    const int nB = self->nB;
    const int* const rise = self->local.rise;
    const int* const last = self->local.last;
    Py_ssize_t* const M_counts = counts;
    Py_ssize_t* const Ix_counts = M_counts + nB + 1;
    Py_ssize_t* const Iy_counts = Ix_counts + nB + 1;
    int j;
    int trace;
    Py_ssize_t count;
    Py_ssize_t term;
    Py_ssize_t M_temp;
    Py_ssize_t Ix_temp;
    Py_ssize_t Iy_temp;
    Py_ssize_t sum = *total;
    M_temp = M_counts[0];
    M_counts[0] = 1;
    Ix_temp = Ix_counts[0];
    Ix_counts[0] = 0;
    Iy_temp = Iy_counts[0];
    Iy_counts[0] = 0;
    for (j = 1; j <= nB; j++) {
        count = 0;
        trace = tr[2*(size_t)j];
        if (trace & ENDPOINT && rise[0] && !linear_before(i, j, rise)
                             && linear_before(i, j, last))
            trace &= ~ENDPOINT;
        if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
        if (trace & Ix_MATRIX) SAFE_ADD(Ix_temp, count);
        if (trace & Iy_MATRIX) SAFE_ADD(Iy_temp, count);
        if (count == 0 && (trace & STARTPOINT)) count = 1;
        M_temp = M_counts[j];
        if (trace & ENDPOINT) {
            SAFE_ADD(count, sum);
            count = 0;
        }
        M_counts[j] = count;
        count = 0;
        trace = LINEAR_GOTOH_IX(tr, j);
        if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
        if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j], count);
        if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j], count);
        Ix_temp = Ix_counts[j];
        Ix_counts[j] = count;
        count = 0;
        trace = LINEAR_GOTOH_IY(tr, j);
        if (trace & M_MATRIX) SAFE_ADD(M_counts[j-1], count);
        if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j-1], count);
        if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j-1], count);
        Iy_temp = Iy_counts[j];
        Iy_counts[j] = count;
    }
    *total = sum;
}

/* The rows of GOTOH_LOCAL_ALIGN: Ix and Iy 0 without a trace in the last
 * row and the last column, and column 0 as there, with the traces
 * PathGenerator_create_Gotoh gives it.  A saved row holds a LinearLocal,
 * then the scores. */
#define LINEAR_GOTOH_LOCAL_SWEEP(align_score) \
    const Aligner* const self = &lp->aligner; \
    const int nA = lp->nA; \
    const int nB = lp->nB; \
    const int* const sA = lp->sA; \
    const int* const sB = lp->sB; \
    const double gap_open_A = self->open_internal_insertion_score; \
    const double gap_open_B = self->open_internal_deletion_score; \
    const double gap_extend_A = self->extend_internal_insertion_score; \
    const double gap_extend_B = self->extend_internal_deletion_score; \
    const double epsilon = self->epsilon; \
    double* const row = saved + LINEAR_LOCAL_DOUBLES; \
    LinearLocal local; \
    double maximum; \
    int i; \
    int j; \
    int jmax; \
    int kA; \
    int kB; \
    int trace; \
    int interrupted = 0; \
    double score; \
    double temp; \
    double M_temp; \
    double Ix_temp; \
    double Iy_temp; \
    unsigned char* tr; \
    memcpy(&local, saved, sizeof(LinearLocal)); \
    maximum = local.maximum; \
    PAIRWISE_NOGIL_BEGIN \
    for (i = r0 + 1; i <= r1; i++) { \
        PAIRWISE_NOGIL_CHECK(i) \
        tr = T + (size_t)(i - r0 - 1) * tstride; \
        M_temp = row[0]; \
        Ix_temp = row[1]; \
        Iy_temp = row[2]; \
        row[0] = 0; \
        row[1] = -DBL_MAX; \
        row[2] = -DBL_MAX; \
        tr[0] = (i == nA) ? 0 : STARTPOINT; \
        tr[1] = 0; \
        kA = sA[i-1]; \
        jmax = (i == nA) ? 0 : (c < nB) ? c : nB - 1; \
        for (j = 1; j <= jmax; j++) { \
            kB = sB[j-1]; \
            LINEAR_GOTOH_LOCAL_CELL(align_score); \
        } \
        for ( ; j <= c; j++) { \
            kB = sB[j-1]; \
            LINEAR_GOTOH_LOCAL_EDGE(align_score); \
        } \
        if (counts) \
            linear_gotoh_local_count(lp, tr, i, counts, &local.total); \
    } \
    PAIRWISE_NOGIL_END \
    local.maximum = maximum; \
    memcpy(saved, &local, sizeof(LinearLocal)); \
    return interrupted ? -1 : 0;

static int
linear_gotoh_local_sweep_compare(const LinearPaths* lp, double* saved, int r0,
                                 int r1, int c, unsigned char* T,
                                 size_t tstride, Py_ssize_t* counts)
{
    const double match = lp->aligner.match;
    const double mismatch = lp->aligner.mismatch;
    const int wildcard = lp->aligner.wildcard;
    LINEAR_GOTOH_LOCAL_SWEEP(COMPARE_SCORE);
}

static int
linear_gotoh_local_sweep_matrix(const LinearPaths* lp, double* saved, int r0,
                                int r1, int c, unsigned char* T,
                                size_t tstride, Py_ssize_t* counts)
{
    const Py_ssize_t n = lp->aligner.substitution_matrix.shape[0];
    LINEAR_GOTOH_LOCAL_SWEEP(MATRIX_SCORE_NOGIL);
}

/* Row 0 of GOTOH_LOCAL_ALIGN, nothing kept yet, and the path counts of row
 * 0 as PathGenerator_gotoh_local_length sets them. */
static void
linear_gotoh_local_start(const LinearPaths* lp, double* saved,
                         Py_ssize_t* counts)
{
    const int nB = lp->nB;
    double* const row = saved + LINEAR_LOCAL_DOUBLES;
    const LinearLocal local = {0};
    int j;
    memcpy(saved, &local, sizeof(LinearLocal));
    for (j = 0; j <= nB; j++) {
        row[3*(size_t)j] = 0;
        row[3*(size_t)j+1] = -DBL_MAX;
        row[3*(size_t)j+2] = -DBL_MAX;
    }
    if (counts) {
        Py_ssize_t* const M_counts = counts;
        Py_ssize_t* const Ix_counts = M_counts + nB + 1;
        Py_ssize_t* const Iy_counts = Ix_counts + nB + 1;
        for (j = 0; j <= nB; j++) {
            M_counts[j] = 1;
            Ix_counts[j] = 0;
            Iy_counts[j] = 0;
        }
    }
}

/* GOTOH_LOCAL_ALIGN's clean-up also removes the traces to matrices that no
 * start point reaches, and any end point there.  No cell before the first
 * end point is one (see linear_sw_score), so before it the clean-up only
 * removes traces to such matrices, which only extreme scores or epsilon,
 * such as infinite or NaN ones, make.  A trace leading to one leads on
 * only to others, and on to an Ix or Iy without a trace, or in row 0 or
 * column 0, where the walk stops and asks the full matrix.  Otherwise the
 * walk, from M, follows the trace bits as the sweep computes them, and
 * finds the full matrix's path.
 *
 * GOTOH_LOCAL_ALIGN raises the highest score only at a rise, so a highest
 * score of 0 can still have end points, scoring within epsilon of it. */
static double
linear_gotoh_local_score(LinearPaths* self, const double* saved)
{
    return linear_local_score(self, saved, M_MATRIX);
}

/* Follow the traceback as PathGenerator_next_gotoh_local does: in each
 * cell, take M if the walker's trace allows it, else Ix, else Iy, to a
 * start point of M.  Row 0 is not in any block: there M is a start point.
 * An empty trace, or Ix or Iy in row 0, is a matrix no start point
 * reaches: the full matrix's clean-up has removed the trace leading to it,
 * and its path goes elsewhere. */
static int
linear_gotoh_local_walk(const LinearPaths* Py_UNUSED(self), LinearWalker* w,
                        const unsigned char* T, int r0, size_t width)
{
    while (w->i > r0) {
        const unsigned char* tr = T + (size_t)(w->i - r0 - 1) * width * 2;
        const size_t j = (size_t)w->j;
        int trace;
        if (w->trace & M_MATRIX) {
            trace = tr[2*j];
            if (trace == STARTPOINT) {
                if (linear_walker_add(w) < 0) return -1;
                return LINEAR_ENDED;
            }
            trace &= M_MATRIX | Ix_MATRIX | Iy_MATRIX;
            if (!trace) return LINEAR_STOPPED;
            if (linear_walker_step(w, DIAGONAL) < 0) return -1;
            w->i--;
            w->j--;
        }
        else if (w->trace & Ix_MATRIX) {
            trace = LINEAR_GOTOH_IX(tr, j);
            if (!trace) return LINEAR_STOPPED;
            if (linear_walker_step(w, VERTICAL) < 0) return -1;
            w->i--;
        }
        else if (w->trace & Iy_MATRIX) {
            trace = LINEAR_GOTOH_IY(tr, j);
            if (!trace) return LINEAR_STOPPED;
            if (linear_walker_step(w, HORIZONTAL) < 0) return -1;
            w->j--;
        }
        else {
            LINEAR_INTERNAL_ERROR;
            return -1;
        }
        w->trace = trace;
    }
    if (w->i > 0) return LINEAR_EXITED;
    if (!(w->trace & M_MATRIX)) return LINEAR_STOPPED;
    if (linear_walker_add(w) < 0) return -1;
    return LINEAR_ENDED;
}

static const LinearKernel linear_gotoh_local_compare = {
    .ncarried = LINEAR_LOCAL_DOUBLES,
    .ndoubles = 3,
    .trace_bytes = 2,
    .ncounts = 3,
    .matrix_cell_bytes = sizeof(Trace) + sizeof(TraceGapsGotoh),
    .matrix_row_pointers = 2,
    .start = linear_gotoh_local_start,
    .sweep = linear_gotoh_local_sweep_compare,
    .score = linear_gotoh_local_score,
    .total = linear_local_total,
    .walk = linear_gotoh_local_walk,
    .align = Aligner_gotoh_local_align_compare,
};

static const LinearKernel linear_gotoh_local_matrix = {
    .ncarried = LINEAR_LOCAL_DOUBLES,
    .ndoubles = 3,
    .trace_bytes = 2,
    .ncounts = 3,
    .matrix_cell_bytes = sizeof(Trace) + sizeof(TraceGapsGotoh),
    .matrix_row_pointers = 2,
    .start = linear_gotoh_local_start,
    .sweep = linear_gotoh_local_sweep_matrix,
    .score = linear_gotoh_local_score,
    .total = linear_local_total,
    .walk = linear_gotoh_local_walk,
    .align = Aligner_gotoh_local_align_matrix,
};


/* The kernel for an alignment, or NULL if it has none. */
static const LinearKernel*
linear_kernel(const Aligner* self, Algorithm algorithm)
{
    const bool matrix = (self->substitution_matrix.obj != NULL);
    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            switch (self->mode) {
                case Global:
                    return matrix ? &linear_nw_matrix : &linear_nw_compare;
                case Local:
                    return matrix ? &linear_sw_matrix : &linear_sw_compare;
                default:
                    return NULL;
            }
        case Gotoh:
            switch (self->mode) {
                case Global:
                    return matrix ? &linear_gotoh_matrix : &linear_gotoh_compare;
                case Local:
                    return matrix ? &linear_gotoh_local_matrix
                                  : &linear_gotoh_local_compare;
                default:
                    return NULL;
            }
        default:
            return NULL;
    }
}


/* -------------------- the strip driver ------------------- */

/* Plan from the limits: as many checkpoints as the checkpoint budget holds,
 * evenly spaced, and strips halved while their trace bits are over the
 * block budget and halving saves memory: it costs a row for the new level,
 * and saves the trace bits of half the rows.  Returns false if the sizes
 * overflow. */
static bool
linear_plan(LinearPlan* plan, const LinearKernel* kernel, int nA, int nB)
{
    const size_t width = (size_t)nB + 1;
    const size_t tb = kernel->trace_bytes;
    size_t rowbytes;
    size_t rows;

    if (width > (SIZE_MAX / sizeof(double) - kernel->ncarried)
                / kernel->ndoubles) return false;
    plan->rowsize = kernel->ncarried + kernel->ndoubles * width;
    rowbytes = plan->rowsize * sizeof(double);
    plan->count = linear_checkpoint_bytes / rowbytes;
    if (plan->count < 1) plan->count = 1;
    if (plan->count > (size_t)nA) plan->count = (size_t)nA;
    plan->stride = ((size_t)nA + plan->count - 1) / plan->count;
    plan->count = ((size_t)nA + plan->stride - 1) / plan->stride;
    plan->depth = 0;
    rows = plan->stride;
    while (rows > 1 && rows > linear_block_bytes / tb / width
                    && rows / 2 > rowbytes / tb / width) {
        rows -= rows / 2;
        plan->depth++;
    }
    if (rows > SIZE_MAX / tb / width || plan->depth > SIZE_MAX / rowbytes)
        return false;
    plan->cells = rows * width;
    /* the sequences, the checkpoints, the row being recomputed, a row per
     * level, a block, and a row of scratch trace bits */
    plan->peak = (double)sizeof(int) * ((double)nA + (double)nB)
               + (double)rowbytes * (double)(plan->count + 1 + plan->depth)
               + (double)tb * ((double)plan->cells + (double)width);
    return true;
}

typedef struct {
    LinearPaths* lp;
    double* work;          /* the rows of a block being recomputed */
    unsigned char* block;  /* the trace bits of a block */
    unsigned char* scratch;/* trace bits of a row, when not kept */
    double* arena;         /* one row per level of bisection */
    LinearWalker walker;
} LinearSolver;

/* Follow the path from the walker, in row r1 and a column c, up to row
 * r0, starting from the values of row r0 in inrow (and what the kernel
 * carries in it).  A strip of at most
 * the planned rows and columns fits in the block by the planned depth, so
 * any narrower or shorter one does too. */
static int
linear_solve(LinearSolver* S, int r0, int r1, const double* inrow,
             size_t level)
{
    const LinearPaths* self = S->lp;
    const LinearKernel* kernel = self->kernel;
    const int c = S->walker.j;
    const size_t rows = (size_t)(r1 - r0);
    const size_t width = (size_t)c + 1;
    const size_t used = kernel->ncarried + kernel->ndoubles * width;
    int status;

    if (rows == 1 || rows <= self->plan.cells / width) {
        memcpy(S->work, inrow, used * sizeof(double));
        if (kernel->sweep(self, S->work, r0, r1, c, S->block,
                          width * kernel->trace_bytes, NULL) < 0) return -1;
        return kernel->walk(self, &S->walker, S->block, r0, width);
    }
    else {
        const int mid = r0 + (int)(rows / 2);
        double* row = S->arena + level * self->plan.rowsize;
        memcpy(row, inrow, used * sizeof(double));
        if (kernel->sweep(self, row, r0, mid, c, S->scratch, 0, NULL) < 0)
            return -1;
        status = linear_solve(S, mid, r1, row, level + 1);
        if (status != LINEAR_EXITED) return status;
        /* the bottom half is done with this level's row */
        return linear_solve(S, r0, mid, inrow, level);
    }
}

/* Find the first path, with the full matrix's tie-breaking.  Returns NULL
 * with an exception set on error, and NULL without one if the path stops
 * before the origin. */
static PyObject*
LinearPaths_first_path(LinearPaths* self)
{
    const LinearKernel* kernel = self->kernel;
    const LinearPlan* plan = &self->plan;
    const size_t rowbytes = plan->rowsize * sizeof(double);
    size_t k;
    int status = LINEAR_EXITED;
    PyObject* path = NULL;
    LinearSolver S = {0};

    S.lp = self;
    S.work = PyMem_Malloc(rowbytes);
    S.block = PyMem_Malloc(plan->cells * kernel->trace_bytes);
    S.scratch = PyMem_Malloc(((size_t)self->nB + 1) * kernel->trace_bytes);
    S.arena = PyMem_Malloc(plan->depth ? plan->depth * rowbytes : 1);
    if (!S.work || !S.block || !S.scratch || !S.arena) {
        PyErr_NoMemory();
        goto exit;
    }
    S.walker.i = self->end_i;
    S.walker.j = self->end_j;
    S.walker.trace = self->end;
    /* the strips from the one holding the end of the path up */
    k = ((size_t)S.walker.i + plan->stride - 1) / plan->stride;
    while (k-- > 0) {
        const size_t r0 = k * plan->stride;
        if ((size_t)S.walker.i <= r0
         || (size_t)S.walker.i > r0 + plan->stride) {
            LINEAR_INTERNAL_ERROR;
            goto exit;
        }
        status = linear_solve(&S, (int)r0, S.walker.i,
                              self->checkpoints + k * plan->rowsize, 0);
        if (status < 0) goto exit;
        if (status != LINEAR_EXITED) break;
    }
    if (status == LINEAR_STOPPED) goto exit;
    if (status != LINEAR_ENDED) {
        LINEAR_INTERNAL_ERROR;
        goto exit;
    }
    path = linear_walker_path(self, &S.walker);
exit:
    PyMem_Free(S.work);
    PyMem_Free(S.block);
    PyMem_Free(S.scratch);
    PyMem_Free(S.arena);
    PyMem_Free(S.walker.points);
    return path;
}

/* The forward pass: compute the score, and save the checkpoints. */
static int
LinearPaths_forward(LinearPaths* self)
{
    const LinearKernel* kernel = self->kernel;
    const LinearPlan* plan = &self->plan;
    const size_t rowbytes = plan->rowsize * sizeof(double);
    size_t k;
    double* row = NULL;
    unsigned char* scratch = NULL;
    int status = -1;

    self->checkpoints = PyMem_Malloc(plan->count * rowbytes);
    row = PyMem_Malloc(rowbytes);
    scratch = PyMem_Malloc(((size_t)self->nB + 1) * kernel->trace_bytes);
    if (!self->checkpoints || !row || !scratch) {
        PyErr_NoMemory();
        goto exit;
    }
    kernel->start(self, row, NULL);
    for (k = 0; k < plan->count; k++) {
        const size_t r0 = k * plan->stride;
        const size_t r1 = Py_MIN(r0 + plan->stride, (size_t)self->nA);
        memcpy(self->checkpoints + k * plan->rowsize, row, rowbytes);
        if (kernel->sweep(self, row, (int)r0, (int)r1, self->nB,
                          scratch, 0, NULL) < 0) goto exit;
    }
    self->score = kernel->score(self, row);
    status = 0;
exit:
    PyMem_Free(row);
    PyMem_Free(scratch);
    return status;
}

/* Count the paths as PathGenerator_length does over the full matrix, in
 * another forward pass over every row, with a row of path counts.  Returns
 * the count, or -1 with an exception set.  The count is kept, and so is an
 * overflow.  A count of 0 is not kept, as 0 means not counted yet, but with
 * no end (see LinearPaths) every kernel's total is 0 without a pass. */
static Py_ssize_t
LinearPaths_count(LinearPaths* self)
{
    const LinearKernel* kernel = self->kernel;
    const size_t width = (size_t)self->nB + 1;
    double* row = NULL;
    unsigned char* scratch = NULL;
    Py_ssize_t* counts = NULL;
    int status = 0;

    if (self->length == 0 && self->end) {
        if (width <= PY_SSIZE_T_MAX / sizeof(Py_ssize_t) / kernel->ncounts) {
            row = PyMem_Malloc(self->plan.rowsize * sizeof(double));
            scratch = PyMem_Malloc(width * kernel->trace_bytes);
            counts = PyMem_Malloc(width * kernel->ncounts * sizeof(Py_ssize_t));
        }
        if (!row || !scratch || !counts) {
            PyErr_NoMemory();
            status = -1;
        }
        else {
            kernel->start(self, row, counts);
            status = kernel->sweep(self, row, 0, self->nA, self->nB,
                                   scratch, 0, counts);
            if (status == 0) self->length = kernel->total(self, row, counts);
        }
        PyMem_Free(row);
        PyMem_Free(scratch);
        PyMem_Free(counts);
        if (status < 0) return -1;
    }
    if (self->length == OVERFLOW_ERROR) {
        PyErr_Format(PyExc_OverflowError,
                     "number of optimal alignments is larger than %zd",
                     PY_SSIZE_T_MAX);
        return -1;
    }
    return self->length;
}

/* Build the full traceback matrix, check it against what is already
 * known, and step its PathGenerator past the first path, which becomes
 * the first path if that is not known yet.  If that is not known, the full
 * matrix of a Gotoh local alignment may have no paths, its clean-up having
 * removed them (see linear_gotoh_local_score); then there is no end.  For
 * any other kernel that is an internal error. */
static int
LinearPaths_materialize(LinearPaths* self)
{
    PyObject* result;
    PyObject* first;
    double score;
    int equal = 1;

    result = self->kernel->align(&self->aligner, self->sA, self->nA,
                                 self->sB, self->nB, self->strand);
    if (!result) return -1;
    score = PyFloat_AsDouble(PyTuple_GET_ITEM(result, 0));
    if (!(score == self->score
          || (Py_IS_NAN(score) && Py_IS_NAN(self->score)))) {
        Py_DECREF(result);
        LINEAR_INTERNAL_ERROR;
        return -1;
    }
    self->paths = (PathGenerator*)PyTuple_GET_ITEM(result, 1);
    Py_INCREF(self->paths);
    Py_DECREF(result);
    first = PathGenerator_next(self->paths);
    if (!first && !self->first && !PyErr_Occurred()
        && (self->kernel == &linear_gotoh_local_compare
            || self->kernel == &linear_gotoh_local_matrix)) {
        self->end = 0;
        PyMem_Free(self->checkpoints);
        self->checkpoints = NULL;
        self->primed = true;
        return 0;
    }
    if (first && self->first)
        equal = PyObject_RichCompareBool(first, self->first, Py_EQ);
    if (!first || equal != 1) {
        if (equal == 0 || (!first && !PyErr_Occurred())) LINEAR_INTERNAL_ERROR;
        Py_XDECREF(first);
        Py_CLEAR(self->paths);
        return -1;
    }
    if (self->first) Py_DECREF(first);
    else {
        self->first = first;
        PyMem_Free(self->checkpoints);
        self->checkpoints = NULL;
    }
    self->primed = true;
    return 0;
}


/* ------------------- the LinearPaths type ------------------- */

static PyObject*
LinearPaths_next(LinearPaths* self)
{
    PyObject* path = NULL;
    if (linear_lock(self) < 0) return NULL;
    if (self->position == 0) {
        if (!self->first) {
            /* no paths, as the full matrix's generator finds at once */
            if (!self->end) goto exit;
            self->first = LinearPaths_first_path(self);
            if (self->first) {
                PyMem_Free(self->checkpoints);
                self->checkpoints = NULL;
            }
            /* a path stopping before the origin, or before a start point,
             * is whatever the full matrix makes of it, if anything */
            else if (PyErr_Occurred() || LinearPaths_materialize(self) < 0
                     || !self->first)
                goto exit;
        }
        Py_INCREF(self->first);
        path = self->first;
    }
    else {
        if (!self->paths && LinearPaths_materialize(self) < 0) goto exit;
        if (!self->primed) {
            /* skip the first path, which was returned already */
            PyObject* first = PathGenerator_next(self->paths);
            if (!first) {
                if (!PyErr_Occurred()) LINEAR_INTERNAL_ERROR;
                goto exit;
            }
            Py_DECREF(first);
            self->primed = true;
        }
        path = PathGenerator_next(self->paths);
        if (!path) goto exit;
    }
    self->position++;
exit:
    linear_unlock(self);
    return path;
}

/* len() counts the paths without building the full matrix, but counts
 * over it once it has been built, as that is quicker. */
static Py_ssize_t
LinearPaths_length(LinearPaths* self)
{
    Py_ssize_t length;
    if (linear_lock(self) < 0) return -1;
    if (self->paths && self->length == 0)
        length = PathGenerator_length(self->paths);
    else length = LinearPaths_count(self);
    linear_unlock(self);
    return length;
}

static PyObject*
LinearPaths_reset(LinearPaths* self, PyObject* Py_UNUSED(ignored))
{
    if (linear_lock(self) < 0) return NULL;
    self->position = 0;
    if (self->paths) {
        Py_DECREF(PathGenerator_reset(self->paths));  /* returns None */
        self->primed = false;
    }
    linear_unlock(self);
    Py_RETURN_NONE;
}

static PyObject*
LinearPaths_get_materialized(LinearPaths* self, void* closure)
{
    return PyBool_FromLong(self->paths != NULL);
}

static void
LinearPaths_dealloc(LinearPaths* self)
{
    PyMem_Free(self->aligner.substitution_matrix.buf);
    PyMem_Free(self->sA);
    PyMem_Free(self->sB);
    PyMem_Free(self->checkpoints);
    Py_XDECREF(self->first);
    Py_XDECREF(self->paths);
    if (self->lock) PyThread_free_lock(self->lock);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMethodDef LinearPaths_methods[] = {
    {"reset", (PyCFunction)LinearPaths_reset, METH_NOARGS,
     PathGenerator_reset__doc__},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static PyGetSetDef LinearPaths_getset[] = {
    {"_materialized", (getter)LinearPaths_get_materialized, NULL,
     "whether the full traceback matrix has been built (for testing)", NULL},
    {NULL}  /* Sentinel */
};

static PySequenceMethods LinearPaths_as_sequence = {
    .sq_length = (lenfunc)LinearPaths_length,
};

static PyTypeObject LinearPaths_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_pairwisealigner.LinearPaths",
    .tp_basicsize = sizeof(LinearPaths),
    .tp_dealloc = (destructor)LinearPaths_dealloc,
    .tp_as_sequence = &LinearPaths_as_sequence,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_iter = PyObject_SelfIter,
    .tp_iternext = (iternextfunc)LinearPaths_next,
    .tp_methods = LinearPaths_methods,
    .tp_getset = LinearPaths_getset,
};


/* ------------------- the hooks into align() ------------------- */

/* Whether align() should return a LinearPaths object: when the full
 * traceback matrix would be over the threshold, and the linear-space
 * traceback would hold less memory than it.  self is the snapshot of the
 * aligner taken by align(). */
static bool
LinearPaths_wanted(const Aligner* self, Algorithm algorithm, int nA, int nB)
{
    const LinearKernel* kernel = linear_kernel(self, algorithm);
    const size_t na = (size_t)nA + 1;
    size_t rowbytes;
    size_t nbytes;  /* of the full traceback matrix */
    double peak;
    LinearPlan plan;
    if (!linear_enabled || !kernel) return false;
    if (!linear_plan(&plan, kernel, nA, nB)) return false;
    if (linear_threshold < 0) return true;
    rowbytes = ((size_t)nB + 1) * kernel->matrix_cell_bytes
             + kernel->matrix_row_pointers * sizeof(void*);
    nbytes = (na <= SIZE_MAX / rowbytes) ? na * rowbytes : SIZE_MAX;
    /* LinearPaths_align keeps a copy of the substitution matrix, on top of
     * the one each kernel makes while it runs without the GIL */
    peak = plan.peak;
    if (self->substitution_matrix.obj)
        peak += (double)self->substitution_matrix.len;
    return nbytes > (size_t)linear_threshold && peak < (double)nbytes;
}

/* Return (score, paths) as the alignment kernels do, for a LinearPaths
 * object.  self is the snapshot of the aligner taken by align(), and the
 * sequences have been mapped to indices already. */
static PyObject*
LinearPaths_align(const Aligner* self, const int* sA, int nA,
                  const int* sB, int nB, unsigned char strand)
{
    /* align() resolved self->algorithm, calling _get_algorithm on its
     * snapshot, before LinearPaths_wanted */
    const LinearKernel* kernel = linear_kernel(self, self->algorithm);
    Py_buffer* view;
    LinearPaths* paths;

    if (!kernel) {
        LINEAR_INTERNAL_ERROR;
        return NULL;
    }
    paths = (LinearPaths*)PyType_GenericAlloc(&LinearPaths_Type, 0);
    if (!paths) return NULL;
    /* PyType_GenericAlloc zeroed every field. */
    paths->kernel = kernel;
    paths->aligner = *self;
    paths->aligner.insertion_score_function = NULL;
    paths->aligner.deletion_score_function = NULL;
    paths->aligner.alphabet = NULL;
    view = &paths->aligner.substitution_matrix;
    memset(view, 0, sizeof(Py_buffer));
    if (self->substitution_matrix.obj) {
        const Py_ssize_t nbytes = self->substitution_matrix.len;
        view->buf = PyMem_Malloc((size_t)nbytes);
        if (!view->buf) goto nomemory;
        memcpy(view->buf, self->substitution_matrix.buf, (size_t)nbytes);
        view->len = nbytes;
        view->itemsize = sizeof(double);
        view->ndim = 2;
        paths->shape[0] = self->substitution_matrix.shape[0];
        paths->shape[1] = self->substitution_matrix.shape[1];
        view->shape = paths->shape;
    }
    paths->nA = nA;
    paths->nB = nB;
    paths->strand = strand;
    if (!linear_plan(&paths->plan, paths->kernel, nA, nB)) goto nomemory;
    paths->sA = PyMem_Malloc((size_t)nA * sizeof(int));
    paths->sB = PyMem_Malloc((size_t)nB * sizeof(int));
    paths->lock = PyThread_allocate_lock();
    if (!paths->sA || !paths->sB || !paths->lock) goto nomemory;
    memcpy(paths->sA, sA, (size_t)nA * sizeof(int));
    memcpy(paths->sB, sB, (size_t)nB * sizeof(int));
    if (LinearPaths_forward(paths) < 0) {
        Py_DECREF(paths);
        return NULL;
    }
    return Py_BuildValue("fN", paths->score, paths);
nomemory:
    Py_DECREF(paths);
    return PyErr_NoMemory();
}

/* _pairwisealigner._set_traceback_limits, for testing. */
static PyObject*
linear_set_traceback_limits(PyObject* module, PyObject* args)
{
    PyObject* threshold;
    Py_ssize_t checkpoint_bytes;
    Py_ssize_t block_bytes;
    Py_ssize_t value = 0;
    PyObject* previous;

    if (!PyArg_ParseTuple(args, "Onn:_set_traceback_limits",
                          &threshold, &checkpoint_bytes, &block_bytes))
        return NULL;
    if (threshold != Py_None) {
        value = PyLong_AsSsize_t(threshold);
        if (value == -1 && PyErr_Occurred()) return NULL;
    }
    if (checkpoint_bytes < 1 || block_bytes < 1) {
        PyErr_SetString(PyExc_ValueError, "the budgets must be at least 1");
        return NULL;
    }
    previous = Py_BuildValue("(Nnn)",
                             linear_enabled ? PyLong_FromSsize_t(linear_threshold)
                                            : (Py_INCREF(Py_None), Py_None),
                             (Py_ssize_t)linear_checkpoint_bytes,
                             (Py_ssize_t)linear_block_bytes);
    if (!previous) return NULL;
    linear_enabled = (threshold != Py_None);
    linear_threshold = value;
    linear_checkpoint_bytes = (size_t)checkpoint_bytes;
    linear_block_bytes = (size_t)block_bytes;
    return previous;
}

static PyMethodDef linear_module_methods[] = {
    {"_set_traceback_limits",
     (PyCFunction)linear_set_traceback_limits,
     METH_VARARGS,
     "_set_traceback_limits(threshold_bytes, checkpoint_bytes, block_bytes)\n"
     "--\n\n"
     "Use the linear-space traceback above threshold_bytes of traceback\n"
     "matrix where it needs less memory (None: never; negative: always);\n"
     "return the previous limits.  Testing only."
    },
    {NULL, NULL, 0, NULL}  /* Sentinel */
};
