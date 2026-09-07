/* Wrapper that compiles the unmodified GSL 1.16 lmder/lmsder sources
 * (lmder.c + lmutil.c + lmpar.c + lmset.c + lmiterate.c + qrsolv.c, covar.c,
 * convergence.c) against a modern GSL, with the 1.16-specific declarations
 * provided here.  GSL 1.16 is Copyright (C) 1996-2013 Brian Gough et al.,
 * GPL-3.0 (see COPYING). */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* Block the GSL 2.x headers whose declarations conflict with the 1.16 code. */
#define __GSL_MULTIFIT_NLIN_H__ 1
#define __GSL_LINALG_H__ 1

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>

#include "g116_lm.h"
#include "g116_linalg.h"
#include "g116_blas.h"   /* inline level-1 BLAS (bit-identical to gslcblas, no DLL calls) */

/* --- 1.16 declarations ---------------------------------------------------- */
#define gsl_multifit_function_fdf g116_function_fdf

typedef struct {
    const char *name;
    size_t size;
    int (*alloc) (void *state, size_t n, size_t p);
    int (*set) (void *state, gsl_multifit_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx);
    int (*iterate) (void *state, gsl_multifit_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx);
    void (*free) (void *state);
} gsl_multifit_fdfsolver_type;

#define GSL_MULTIFIT_FN_EVAL_F(F,x,y) ((*((F)->f))(x,(F)->params,(y)))
#define GSL_MULTIFIT_FN_EVAL_DF(F,x,dy) ((*((F)->df))(x,(F)->params,(dy)))
#define GSL_MULTIFIT_FN_EVAL_F_DF(F,x,y,dy) ((*((F)->fdf))(x,(F)->params,(y),(dy)))

/* numerical-derivative fallbacks are never used (df/fdf are always provided) */
static int g116_dif_df_stub(const gsl_vector * x, gsl_multifit_function_fdf * fdf, const gsl_vector * f, gsl_matrix * J)
{ (void)x; (void)fdf; (void)f; (void)J; return GSL_EUNIMPL; }
static int g116_dif_fdf_stub(const gsl_vector * x, gsl_multifit_function_fdf * fdf, gsl_vector * f, gsl_matrix * J)
{ (void)x; (void)fdf; (void)f; (void)J; return GSL_EUNIMPL; }
#define gsl_multifit_fdfsolver_dif_df g116_dif_df_stub
#define gsl_multifit_fdfsolver_dif_fdf g116_dif_fdf_stub

/* optional phase timing (cme_detect --profile); the wrappers only add a
 * branch per call when profiling is off */
int cme_prof_enabled_c(void);
double cme_prof_now_c(void);
void cme_prof_add_c(int phase, double t0);
#define CME_PROF_LM_QRPT 19    /* cme::prof::LM_qrpt  (checked by static_assert in profile.cpp) */
#define CME_PROF_LM_QTVEC 20   /* cme::prof::LM_qtvec */
static int g116_QRPT_decomp_timed(gsl_matrix * A, gsl_vector * tau, gsl_permutation * p, int *signum, gsl_vector * norm)
{
    if (cme_prof_enabled_c()) {
        double t0 = cme_prof_now_c();
        int st = g116_QRPT_decomp(A, tau, p, signum, norm);
        cme_prof_add_c(CME_PROF_LM_QRPT, t0);
        return st;
    }
    return g116_QRPT_decomp(A, tau, p, signum, norm);
}
static int g116_QR_QTvec_timed(const gsl_matrix * QR, const gsl_vector * tau, gsl_vector * v)
{
    if (cme_prof_enabled_c()) {
        double t0 = cme_prof_now_c();
        int st = g116_QR_QTvec(QR, tau, v);
        cme_prof_add_c(CME_PROF_LM_QTVEC, t0);
        return st;
    }
    return g116_QR_QTvec(QR, tau, v);
}

/* renames */
#define gsl_linalg_QRPT_decomp g116_QRPT_decomp_timed
#define gsl_linalg_QR_QTvec g116_QR_QTvec_timed
#define gsl_multifit_covar g116_covar
#define gsl_multifit_test_delta g116_test_delta
#define gsl_multifit_test_gradient g116_test_gradient
#define gsl_multifit_fdfsolver_lmder g116_type_lmder
#define gsl_multifit_fdfsolver_lmsder g116_type_lmsder
#define gsl_multifit_gradient g116_gradient_unused

/* config.h is included by the 1.16 sources; provide an empty one via the
 * include path (third_party/gsl116/config.h). */
#include "lmder.c"
#include "covar.c"
#include "convergence.c"

/* --- driver ----------------------------------------------------------------- */
struct g116_solver {
    const gsl_multifit_fdfsolver_type *type;
    void *state;
    gsl_vector *x;
    gsl_vector *f;
    gsl_matrix *J;
    gsl_vector *dx;
};

g116_solver *g116_solver_alloc(size_t n, size_t p)
{
    g116_solver *s = (g116_solver *) calloc(1, sizeof(g116_solver));
    if (!s) return NULL;
    s->type = &lmsder_type;
    s->x = gsl_vector_calloc(p);
    s->f = gsl_vector_calloc(n);
    s->J = gsl_matrix_calloc(n, p);
    s->dx = gsl_vector_calloc(p);
    s->state = calloc(1, s->type->size);
    if (!s->x || !s->f || !s->J || !s->dx || !s->state) { g116_solver_free(s); return NULL; }
    if (s->type->alloc(s->state, n, p) != GSL_SUCCESS) { g116_solver_free(s); return NULL; }
    return s;
}

int g116_solver_set(g116_solver * s, g116_function_fdf * fdf, const gsl_vector * x0)
{
    gsl_vector_memcpy(s->x, x0);
    return s->type->set(s->state, fdf, s->x, s->f, s->J, s->dx);
}

int g116_solver_iterate(g116_solver * s, g116_function_fdf * fdf)
{
    return s->type->iterate(s->state, fdf, s->x, s->f, s->J, s->dx);
}

const gsl_vector *g116_solver_x(const g116_solver * s) { return s->x; }
const gsl_vector *g116_solver_f(const g116_solver * s) { return s->f; }
const gsl_vector *g116_solver_dx(const g116_solver * s) { return s->dx; }
const gsl_matrix *g116_solver_J(const g116_solver * s) { return s->J; }

void g116_solver_free(g116_solver * s)
{
    if (!s) return;
    if (s->state) { s->type->free(s->state); free(s->state); }
    if (s->x) gsl_vector_free(s->x);
    if (s->f) gsl_vector_free(s->f);
    if (s->J) gsl_matrix_free(s->J);
    if (s->dx) gsl_vector_free(s->dx);
    free(s);
}
