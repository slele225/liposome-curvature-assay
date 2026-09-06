/* Public C API of the GSL 1.16 scaled Levenberg-Marquardt solver (lmsder)
 * compiled from the original GSL 1.16 sources (GPL-3.0, see COPYING).
 *
 * The fitGaussian2D MEX of cmeAnalysis was built against GSL 1.x; the
 * lmsder implementation changed in GSL 2.x (column scaling via dnrm2, and the
 * 1.x update_diag() quirk that only sums the first p rows of J), which alters
 * the iteration path of ill-conditioned fits.  Using the 1.16 code makes the
 * port follow the MEX iteration path exactly. */
#ifndef G116_LM_H
#define G116_LM_H
#include <stddef.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int (*f)(const gsl_vector * x, void *params, gsl_vector * f);
    int (*df)(const gsl_vector * x, void *params, gsl_matrix * df);
    int (*fdf)(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * df);
    size_t n;   /* number of functions */
    size_t p;   /* number of independent variables */
    void *params;
} g116_function_fdf;

typedef struct g116_solver g116_solver;

g116_solver *g116_solver_alloc(size_t n, size_t p);           /* lmsder */
int g116_solver_set(g116_solver * s, g116_function_fdf * fdf, const gsl_vector * x0);
int g116_solver_iterate(g116_solver * s, g116_function_fdf * fdf);
const gsl_vector *g116_solver_x(const g116_solver * s);
const gsl_vector *g116_solver_f(const g116_solver * s);
const gsl_vector *g116_solver_dx(const g116_solver * s);
const gsl_matrix *g116_solver_J(const g116_solver * s);
void g116_solver_free(g116_solver * s);

int g116_test_delta(const gsl_vector * dx, const gsl_vector * x, double epsabs, double epsrel);
int g116_covar(const gsl_matrix * J, double epsrel, gsl_matrix * covar);

#ifdef __cplusplus
}
#endif
#endif
