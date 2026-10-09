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
 * - The first path is recovered strip by strip from the bottom: recompute
 *   a strip's rows from its checkpoint, keeping their trace bits, and
 *   follow the traceback up through it.  Only the columns left of where the
 *   path leaves the strip are needed.  A strip over the block budget is
 *   bisected where that saves memory: sweep to its middle row, solve the
 *   bottom half, then the top half.  Same per-cell code, same trace bits,
 *   same path.  (Hirschberg's forward-and-reverse midpoints would pick
 *   another co-optimal path.)
 * - len() counts the paths in one more forward pass, recomputing each row's
 *   trace bits and keeping a row of counts.
 * - Later paths build the full matrix with the existing kernel and delegate
 *   to its PathGenerator, after checking it agrees.
 *
 * It is used only where it holds less memory than the matrix.  A row of
 * trace bits costs less than a row of doubles, so a short target against
 * a long query has too few rows to save.
 *
 * Settings, substitution matrix and sequences are copied at align(), so
 * later changes to them cannot change the results.  Each algorithm
 * supplies a LinearKernel to the generic strip driver; only
 * Needleman-Wunsch has one so far.  The threshold is 512 MiB by default;
 * the private _set_traceback_limits changes it, for testing.
 */


typedef struct LinearPaths LinearPaths;
typedef struct LinearWalker LinearWalker;

/* Return values of LinearKernel.walk and of linear_solve. */
#define LINEAR_EXITED 0   /* the walk reached the top row of the block */
#define LINEAR_ENDED 1    /* the path is complete */

typedef struct {
    /* A saved row holds ncarried scalars, then ndoubles per column. */
    size_t ncarried;
    size_t ndoubles;
    size_t trace_bytes;  /* bytes of trace bits per cell */
    size_t ncounts;      /* path counts per column, for len() */
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
    /* The score of the alignment, from the last row. */
    double (*score)(const LinearPaths* self, const double* row);
    /* The number of paths, from the path counts of the last row, or
     * OVERFLOW_ERROR. */
    Py_ssize_t (*total)(const LinearPaths* self, const Py_ssize_t* counts);
    /* Follow the traceback up through the block holding the trace bits of
     * rows r0+1 to the walker's row, width cells per row.  Returns
     * LINEAR_EXITED with the walker in row r0, LINEAR_ENDED, or -1 with an
     * exception set. */
    int (*walk)(const LinearPaths* self, LinearWalker* walker,
                const unsigned char* T, int r0, size_t width);
    /* The kernel building the full traceback matrix. */
    PyObject* (*align)(Aligner* self, const int* sA, int nA,
                       const int* sB, int nB, unsigned char strand);
} LinearKernel;

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
linear_nw_score(const LinearPaths* self, const double* row)
{
    return row[self->nB];
}

static Py_ssize_t
linear_nw_total(const LinearPaths* self, const Py_ssize_t* counts)
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
    .start = linear_nw_start,
    .sweep = linear_nw_sweep_matrix,
    .score = linear_nw_score,
    .total = linear_nw_total,
    .walk = linear_nw_walk,
    .align = Aligner_needlemanwunsch_align_matrix,
};


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
 * r0, starting from the values of row r0 in inrow.  A strip of at most
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

/* Find the first path, with the full matrix's tie-breaking. */
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
    S.walker.i = self->nA;
    S.walker.j = self->nB;
    for (k = plan->count; k-- > 0; ) {
        const size_t r0 = k * plan->stride;
        const size_t r1 = Py_MIN(r0 + plan->stride, (size_t)self->nA);
        if (S.walker.i != (int)r1) {
            LINEAR_INTERNAL_ERROR;
            goto exit;
        }
        status = linear_solve(&S, (int)r0, (int)r1,
                              self->checkpoints + k * plan->rowsize, 0);
        if (status < 0) goto exit;
        if (status == LINEAR_ENDED) break;
    }
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
 * overflow. */
static Py_ssize_t
LinearPaths_count(LinearPaths* self)
{
    const LinearKernel* kernel = self->kernel;
    const size_t width = (size_t)self->nB + 1;
    double* row = NULL;
    unsigned char* scratch = NULL;
    Py_ssize_t* counts = NULL;
    int status = 0;

    if (self->length == 0) {
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
            if (status == 0) self->length = kernel->total(self, counts);
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
 * known, and step its PathGenerator past the first path. */
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
            self->first = LinearPaths_first_path(self);
            if (!self->first) goto exit;
            PyMem_Free(self->checkpoints);
            self->checkpoints = NULL;
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
    const size_t na = (size_t)nA + 1;
    const size_t rowbytes = ((size_t)nB + 1) * sizeof(Trace) + sizeof(Trace*);
    size_t nbytes;  /* of the full traceback matrix */
    double peak;
    LinearPlan plan;
    if (!linear_enabled) return false;
    if (algorithm != NeedlemanWunschSmithWaterman || self->mode != Global)
        return false;
    /* both Needleman-Wunsch kernels have the same sizes */
    if (!linear_plan(&plan, &linear_nw_compare, nA, nB)) return false;
    if (linear_threshold < 0) return true;
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
    Py_buffer* view;
    LinearPaths* paths;

    paths = (LinearPaths*)PyType_GenericAlloc(&LinearPaths_Type, 0);
    if (!paths) return NULL;
    /* PyType_GenericAlloc zeroed every field. */
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
        paths->kernel = &linear_nw_matrix;
    }
    else paths->kernel = &linear_nw_compare;
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
