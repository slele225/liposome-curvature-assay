/* Prototypes for the GSL 1.16 linear-algebra extracts (g116_linalg.c). */
#ifndef G116_LINALG_H
#define G116_LINALG_H
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

#ifdef __cplusplus
extern "C" {
#endif

double g116_householder_transform (gsl_vector * v);
int g116_householder_hm (double tau, const gsl_vector * v, gsl_matrix * A);
int g116_householder_hv (double tau, const gsl_vector * v, gsl_vector * w);
int g116_QRPT_decomp (gsl_matrix * A, gsl_vector * tau, gsl_permutation * p, int *signum, gsl_vector * norm);
int g116_QR_QTvec (const gsl_matrix * QR, const gsl_vector * tau, gsl_vector * v);

#ifdef __cplusplus
}
#endif
#endif
