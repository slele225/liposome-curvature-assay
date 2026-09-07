/* Inline replacements for the handful of level-1 BLAS wrappers used by the
 * GSL 1.16 Levenberg-Marquardt code (lmder/lmpar/qrsolv/householder).
 *
 * The vendored GSL is a DLL, so every gsl_blas_* call crosses two DLL
 * boundaries (gsl.dll -> gslcblas.dll) for vectors of a few hundred elements.
 * These functions reproduce the reference gslcblas implementations
 * (cblas/source_dot_r.h, source_axpy_r.h, source_scal_r.h, source_nrm2_r.h,
 * source_iamax_r.h) operation by operation, in the same order, so the
 * results are bit-identical; only the call overhead is removed.  The
 * accumulation loops are plain sequential loops (no reassociation; the code
 * is compiled with /fp:precise, which forbids FP contraction and
 * re-association).
 *
 * Copyright of the reference algorithms: GSL, GPL-3.0 (see COPYING). */
#ifndef G116_BLAS_H
#define G116_BLAS_H
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>

static double g116_cblas_ddot(size_t N, const double *X, size_t incX, const double *Y, size_t incY)
{
    double r = 0.0;
    size_t i, ix = 0, iy = 0;
    for (i = 0; i < N; i++) {
        r += X[ix] * Y[iy];
        ix += incX;
        iy += incY;
    }
    return r;
}

static void g116_cblas_daxpy(size_t N, double alpha, const double *X, size_t incX, double *Y, size_t incY)
{
    size_t i, ix = 0, iy = 0;
    if (alpha == 0.0) return;
    for (i = 0; i < N; i++) {
        Y[iy] += alpha * X[ix];
        ix += incX;
        iy += incY;
    }
}

static void g116_cblas_dscal(size_t N, double alpha, double *X, size_t incX)
{
    size_t i, ix = 0;
    for (i = 0; i < N; i++) {
        X[ix] *= alpha;
        ix += incX;
    }
}

static double g116_cblas_dnrm2(size_t N, const double *X, size_t incX)
{
    double scale = 0.0;
    double ssq = 1.0;
    size_t i, ix = 0;
    if (N == 0) return 0.0;
    if (N == 1) return fabs(X[0]);
    /* the reference writes (scale/ax)*(scale/ax) and (ax/scale)*(ax/scale);
     * the two divisions are identical, so computing the quotient once gives
     * the same bits and halves the division count */
    for (i = 0; i < N; i++) {
        const double x = X[ix];
        if (x != 0.0) {
            const double ax = fabs(x);
            if (scale < ax) {
                const double t = scale / ax;
                ssq = 1.0 + ssq * t * t;
                scale = ax;
            } else {
                const double t = ax / scale;
                ssq += t * t;
            }
        }
        ix += incX;
    }
    return scale * sqrt(ssq);
}

static size_t g116_cblas_idamax(size_t N, const double *X, size_t incX)
{
    double max = 0.0;
    size_t i, ix = 0, result = 0;
    for (i = 0; i < N; i++) {
        if (fabs(X[ix]) > max) {
            max = fabs(X[ix]);
            result = i;
        }
        ix += incX;
    }
    return result;
}

/* gsl_vector wrappers with the same argument checks as blas/blas.c */
static int g116_blas_ddot(const gsl_vector *X, const gsl_vector *Y, double *result)
{
    if (X->size != Y->size) GSL_ERROR("invalid length", GSL_EBADLEN);
    *result = g116_cblas_ddot(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
}
static int g116_blas_daxpy(double alpha, const gsl_vector *X, gsl_vector *Y)
{
    if (X->size != Y->size) GSL_ERROR("invalid length", GSL_EBADLEN);
    g116_cblas_daxpy(X->size, alpha, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
}
static void g116_blas_dscal(double alpha, gsl_vector *X)
{
    g116_cblas_dscal(X->size, alpha, X->data, X->stride);
}
static double g116_blas_dnrm2(const gsl_vector *X)
{
    return g116_cblas_dnrm2(X->size, X->data, X->stride);
}
static size_t g116_blas_idamax(const gsl_vector *X)
{
    return g116_cblas_idamax(X->size, X->data, X->stride);
}

#define gsl_blas_ddot g116_blas_ddot
#define gsl_blas_daxpy g116_blas_daxpy
#define gsl_blas_dscal g116_blas_dscal
#define gsl_blas_dnrm2 g116_blas_dnrm2
#define gsl_blas_idamax g116_blas_idamax

/* Pure data-movement helpers (copies, scaling by a constant, permutation)
 * that the DLL otherwise provides.  No arithmetic is reordered: a copy is a
 * copy, x *= -1 is exact, and applying a permutation moves elements. */
static int g116_vector_memcpy(gsl_vector *dest, const gsl_vector *src)
{
    const size_t n = src->size;
    size_t i;
    if (dest->size != n) GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
    if (dest->stride == 1 && src->stride == 1) {
        for (i = 0; i < n; i++) dest->data[i] = src->data[i];
    } else {
        for (i = 0; i < n; i++) dest->data[i * dest->stride] = src->data[i * src->stride];
    }
    return GSL_SUCCESS;
}
static int g116_matrix_memcpy(gsl_matrix *dest, const gsl_matrix *src)
{
    const size_t n1 = src->size1, n2 = src->size2;
    size_t i, j;
    if (dest->size1 != n1 || dest->size2 != n2) GSL_ERROR("matrix sizes are different", GSL_EBADLEN);
    for (i = 0; i < n1; i++)
        for (j = 0; j < n2; j++)
            dest->data[i * dest->tda + j] = src->data[i * src->tda + j];
    return GSL_SUCCESS;
}
static int g116_vector_scale(gsl_vector *a, const double x)
{
    const size_t n = a->size, stride = a->stride;
    size_t i;
    for (i = 0; i < n; i++) a->data[i * stride] *= x;
    return GSL_SUCCESS;
}
static void g116_vector_set_all(gsl_vector *v, double x)
{
    const size_t n = v->size, stride = v->stride;
    size_t i;
    for (i = 0; i < n; i++) v->data[i * stride] = x;
}
static void g116_vector_set_zero(gsl_vector *v)
{
    g116_vector_set_all(v, 0.0);
}
/* gsl_permute_vector_inverse: v'[p[i]] = v[i] */
static int g116_permute_vector_inverse(const gsl_permutation *p, gsl_vector *v)
{
    const size_t n = v->size, stride = v->stride;
    double tmp_stack[16];
    double *tmp = n <= 16 ? tmp_stack : (double *)malloc(n * sizeof(double));
    size_t i;
    if (p->size != n) GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
    for (i = 0; i < n; i++) tmp[p->data[i]] = v->data[i * stride];
    for (i = 0; i < n; i++) v->data[i * stride] = tmp[i];
    if (tmp != tmp_stack) free(tmp);
    return GSL_SUCCESS;
}
#define gsl_vector_memcpy g116_vector_memcpy
#define gsl_matrix_memcpy g116_matrix_memcpy
#define gsl_vector_scale g116_vector_scale
#define gsl_vector_set_all g116_vector_set_all
#define gsl_vector_set_zero g116_vector_set_zero
#define gsl_permute_vector_inverse g116_permute_vector_inverse

#endif
