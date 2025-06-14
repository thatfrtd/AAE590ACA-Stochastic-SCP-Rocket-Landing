/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes.c
 *
 * Code generation for function 'mtimes'
 *
 */

/* Include files */
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "blas.h"
#include <stddef.h>

/* Function Definitions */
void b_mtimes(const real_T A[3000], const real_T B[3], real_T C[1000])
{
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  real_T alpha1;
  real_T beta1;
  char_T TRANSA1;
  char_T TRANSB1;
  TRANSB1 = 'N';
  TRANSA1 = 'T';
  alpha1 = 1.0;
  beta1 = 0.0;
  m_t = (ptrdiff_t)1000;
  n_t = (ptrdiff_t)1;
  k_t = (ptrdiff_t)3;
  lda_t = (ptrdiff_t)3;
  ldb_t = (ptrdiff_t)3;
  ldc_t = (ptrdiff_t)1000;
  dgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &alpha1, (real_T *)&A[0], &lda_t,
        (real_T *)&B[0], &ldb_t, &beta1, &C[0], &ldc_t);
}

void c_mtimes(const real_T A[3000], const real_T B[1000], real_T C[3])
{
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  real_T alpha1;
  real_T beta1;
  char_T TRANSA1;
  char_T TRANSB1;
  TRANSB1 = 'N';
  TRANSA1 = 'N';
  alpha1 = 1.0;
  beta1 = 0.0;
  m_t = (ptrdiff_t)3;
  n_t = (ptrdiff_t)1;
  k_t = (ptrdiff_t)1000;
  lda_t = (ptrdiff_t)3;
  ldb_t = (ptrdiff_t)1000;
  ldc_t = (ptrdiff_t)3;
  dgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &alpha1, (real_T *)&A[0], &lda_t,
        (real_T *)&B[0], &ldb_t, &beta1, &C[0], &ldc_t);
}

void mtimes(const real_T A[1000000], const real_T B[1000], real_T C[1000])
{
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  real_T alpha1;
  real_T beta1;
  char_T TRANSA1;
  char_T TRANSB1;
  TRANSB1 = 'N';
  TRANSA1 = 'N';
  alpha1 = 1.0;
  beta1 = 0.0;
  m_t = (ptrdiff_t)1000;
  n_t = (ptrdiff_t)1;
  k_t = (ptrdiff_t)1000;
  lda_t = (ptrdiff_t)1000;
  ldb_t = (ptrdiff_t)1000;
  ldc_t = (ptrdiff_t)1000;
  dgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &alpha1, (real_T *)&A[0], &lda_t,
        (real_T *)&B[0], &ldb_t, &beta1, &C[0], &ldc_t);
}

/* End of code generation (mtimes.c) */
