/* Verbatim extracts from GSL 1.16 (linalg/householder.c, linalg/qrpt.c,
 * linalg/qr.c), Copyright (C) 1996-2007 Gerard Jungman, Brian Gough,
 * GPL-3.0 (see COPYING).  Renamed with the g116_ prefix so that the
 * Levenberg-Marquardt solver of GSL 1.16 can be reproduced bit-for-bit next to
 * a modern GSL. */
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include "g116_linalg.h"
#include "g116_blas.h"   /* inline level-1 BLAS (bit-identical to gslcblas, no DLL calls) */

double
g116_householder_transform (gsl_vector * v)
{
  /* replace v[0:n-1] with a householder vector (v[0:n-1]) and
     coefficient tau that annihilate v[1:n-1] */

  const size_t n = v->size ;

  if (n == 1)
    {
      return 0.0; /* tau = 0 */
    }
  else
    { 
      double alpha, beta, tau ;
      
      gsl_vector_view x = gsl_vector_subvector (v, 1, n - 1) ; 
      
      double xnorm = gsl_blas_dnrm2 (&x.vector);
      
      if (xnorm == 0) 
        {
          return 0.0; /* tau = 0 */
        }
      
      alpha = gsl_vector_get (v, 0) ;
      beta = - (alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, xnorm) ;
      tau = (beta - alpha) / beta ;
      
      {
        double s = (alpha - beta);
        
        if (fabs(s) > GSL_DBL_MIN) 
          {
            gsl_blas_dscal (1.0 / s, &x.vector);
            gsl_vector_set (v, 0, beta) ;
          }
        else
          {
            gsl_blas_dscal (GSL_DBL_EPSILON / s, &x.vector);
            gsl_blas_dscal (1.0 / GSL_DBL_EPSILON, &x.vector);
            gsl_vector_set (v, 0, beta) ;
          }
      }
      
      return tau;
    }
}

int
g116_householder_hm (double tau, const gsl_vector * v, gsl_matrix * A)
{
  /* applies a householder transformation v,tau to matrix m */

  if (tau == 0.0)
    {
      return GSL_SUCCESS;
    }

#ifdef USE_BLAS
  {
    gsl_vector_const_view v1 = gsl_vector_const_subvector (v, 1, v->size - 1);
    gsl_matrix_view A1 = gsl_matrix_submatrix (A, 1, 0, A->size1 - 1, A->size2);
    size_t j;

    for (j = 0; j < A->size2; j++)
      {
        double wj = 0.0;
        gsl_vector_view A1j = gsl_matrix_column(&A1.matrix, j);
        gsl_blas_ddot (&A1j.vector, &v1.vector, &wj);
        wj += gsl_matrix_get(A,0,j);

        {
          double A0j = gsl_matrix_get (A, 0, j);
          gsl_matrix_set (A, 0, j, A0j - tau *  wj);
        }

        gsl_blas_daxpy (-tau * wj, &v1.vector, &A1j.vector);
      }
  }
#else
  {
    size_t i, j;
    
    for (j = 0; j < A->size2; j++)
      {
        /* Compute wj = Akj vk */
        
        double wj = gsl_matrix_get(A,0,j);  
        
        for (i = 1; i < A->size1; i++)  /* note, computed for v(0) = 1 above */
          {
            wj += gsl_matrix_get(A,i,j) * gsl_vector_get(v,i);
          }
        
        /* Aij = Aij - tau vi wj */
        
        /* i = 0 */
        {
          double A0j = gsl_matrix_get (A, 0, j);
          gsl_matrix_set (A, 0, j, A0j - tau *  wj);
        }
        
        /* i = 1 .. M-1 */
        
        for (i = 1; i < A->size1; i++)
          {
            double Aij = gsl_matrix_get (A, i, j);
            double vi = gsl_vector_get (v, i);
            gsl_matrix_set (A, i, j, Aij - tau * vi * wj);
          }
      }
  }
#endif
    
  return GSL_SUCCESS;
}

int
g116_householder_hv (double tau, const gsl_vector * v, gsl_vector * w)
{
  /* applies a householder transformation v to vector w */
  const size_t N = v->size;
 
  if (tau == 0)
    return GSL_SUCCESS ;

  {
    /* compute d = v'w */

    double d0 = gsl_vector_get(w,0);
    double d1, d;

    gsl_vector_const_view v1 = gsl_vector_const_subvector(v, 1, N-1);
    gsl_vector_view w1 = gsl_vector_subvector(w, 1, N-1);

    gsl_blas_ddot (&v1.vector, &w1.vector, &d1);
    
    d = d0 + d1;

    /* compute w = w - tau (v) (v'w) */
  
    {
      double w0 = gsl_vector_get (w,0);
      gsl_vector_set (w, 0, w0 - tau * d);
    }
    
    gsl_blas_daxpy (-tau * d, &v1.vector, &w1.vector);
  }
  
  return GSL_SUCCESS;
}

int
g116_QRPT_decomp (gsl_matrix * A, gsl_vector * tau, gsl_permutation * p, int *signum, gsl_vector * norm)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (tau->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (p->size != N)
    {
      GSL_ERROR ("permutation size must be N", GSL_EBADLEN);
    }
  else if (norm->size != N)
    {
      GSL_ERROR ("norm size must be N", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      *signum = 1;

      gsl_permutation_init (p); /* set to identity */

      /* Compute column norms and store in workspace */

      for (i = 0; i < N; i++)
        {
          gsl_vector_view c = gsl_matrix_column (A, i);
          double x = gsl_blas_dnrm2 (&c.vector);
          gsl_vector_set (norm, i, x);
        }

      for (i = 0; i < GSL_MIN (M, N); i++)
        {
          /* Bring the column of largest norm into the pivot position */

          double max_norm = gsl_vector_get(norm, i);
          size_t j, kmax = i;

          for (j = i + 1; j < N; j++)
            {
              double x = gsl_vector_get (norm, j);

              if (x > max_norm)
                {
                  max_norm = x;
                  kmax = j;
                }
            }

          if (kmax != i)
            {
              gsl_matrix_swap_columns (A, i, kmax);
              gsl_permutation_swap (p, i, kmax);
              gsl_vector_swap_elements(norm,i,kmax);

              (*signum) = -(*signum);
            }

          /* Compute the Householder transformation to reduce the j-th
             column of the matrix to a multiple of the j-th unit vector */

          {
            gsl_vector_view c_full = gsl_matrix_column (A, i);
            gsl_vector_view c = gsl_vector_subvector (&c_full.vector, 
                                                      i, M - i);
            double tau_i = g116_householder_transform (&c.vector);

            gsl_vector_set (tau, i, tau_i);

            /* Apply the transformation to the remaining columns */

            if (i + 1 < N)
              {
                gsl_matrix_view m = gsl_matrix_submatrix (A, i, i + 1, M - i, N - (i+1));

                g116_householder_hm (tau_i, &c.vector, &m.matrix);
              }
          }

          /* Update the norms of the remaining columns too */

          if (i + 1 < M) 
            {
              for (j = i + 1; j < N; j++)
                {
                  double x = gsl_vector_get (norm, j);

                  if (x > 0.0)
                    {
                      double y = 0;
                      double temp= gsl_matrix_get (A, i, j) / x;
                  
                      if (fabs (temp) >= 1)
                        y = 0.0;
                      else
                        y = x * sqrt (1 - temp * temp);
                      
                      /* recompute norm to prevent loss of accuracy */

                      if (fabs (y / x) < sqrt (20.0) * GSL_SQRT_DBL_EPSILON)
                        {
                          gsl_vector_view c_full = gsl_matrix_column (A, j);
                          gsl_vector_view c = 
                            gsl_vector_subvector(&c_full.vector,
                                                 i+1, M - (i+1));
                          y = gsl_blas_dnrm2 (&c.vector);
                        }
                  
                      gsl_vector_set (norm, j, y);
                    }
                }
            }
        }

      return GSL_SUCCESS;
    }
}

int
g116_QR_QTvec (const gsl_matrix * QR, const gsl_vector * tau, gsl_vector * v)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (tau->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (v->size != M)
    {
      GSL_ERROR ("vector size must be M", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* compute Q^T v */

      for (i = 0; i < GSL_MIN (M, N); i++)
        {
          gsl_vector_const_view c = gsl_matrix_const_column (QR, i);
          gsl_vector_const_view h = gsl_vector_const_subvector (&(c.vector), i, M - i);
          gsl_vector_view w = gsl_vector_subvector (v, i, M - i);
          double ti = gsl_vector_get (tau, i);
          g116_householder_hv (ti, &(h.vector), &(w.vector));
        }
      return GSL_SUCCESS;
    }
}
