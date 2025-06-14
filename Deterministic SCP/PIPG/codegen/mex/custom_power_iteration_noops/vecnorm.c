/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * vecnorm.c
 *
 * Code generation for function 'vecnorm'
 *
 */

/* Include files */
#include "vecnorm.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
void b_vecnorm(const real_T x[30], real_T y[15])
{
  int32_T j;
  int32_T k;
  for (j = 0; j < 15; j++) {
    real_T b_y;
    real_T scale;
    int32_T ix0;
    int32_T kend;
    ix0 = j << 1;
    b_y = 0.0;
    scale = 3.3121686421112381E-170;
    kend = ix0 + 2;
    for (k = ix0 + 1; k <= kend; k++) {
      real_T absxk;
      absxk = muDoubleScalarAbs(x[k - 1]);
      if (absxk > scale) {
        real_T t;
        t = scale / absxk;
        b_y = b_y * t * t + 1.0;
        scale = absxk;
      } else {
        real_T t;
        t = absxk / scale;
        b_y += t * t;
      }
    }
    y[j] = scale * muDoubleScalarSqrt(b_y);
  }
}

void vecnorm(const real_T x[105], real_T y[15])
{
  int32_T j;
  int32_T k;
  for (j = 0; j < 15; j++) {
    real_T b_y;
    real_T scale;
    int32_T ix0;
    int32_T kend;
    ix0 = j * 7;
    b_y = 0.0;
    scale = 3.3121686421112381E-170;
    kend = ix0 + 7;
    for (k = ix0 + 1; k <= kend; k++) {
      real_T absxk;
      absxk = muDoubleScalarAbs(x[k - 1]);
      if (absxk > scale) {
        real_T t;
        t = scale / absxk;
        b_y = b_y * t * t + 1.0;
        scale = absxk;
      } else {
        real_T t;
        t = absxk / scale;
        b_y += t * t;
      }
    }
    y[j] = scale * muDoubleScalarSqrt(b_y);
  }
}

/* End of code generation (vecnorm.c) */
