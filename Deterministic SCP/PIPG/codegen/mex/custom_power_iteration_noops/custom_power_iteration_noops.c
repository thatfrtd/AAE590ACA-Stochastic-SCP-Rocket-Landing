/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * custom_power_iteration_noops.c
 *
 * Code generation for function 'custom_power_iteration_noops'
 *
 */

/* Include files */
#include "custom_power_iteration_noops.h"
#include "custom_power_iteration_noops_data.h"
#include "rt_nonfinite.h"
#include "sumMatrixIncludeNaN.h"
#include "vecnorm.h"
#include "mwmathutil.h"
#include <emmintrin.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo b_emlrtRSI = {
    27,                             /* lineNo */
    "custom_power_iteration_noops", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\custom_power_iteration_noops.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    52,                             /* lineNo */
    "custom_power_iteration_noops", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\custom_power_iteration_noops.m" /* pathName */
};

static emlrtRTEInfo emlrtRTEI = {
    13,     /* lineNo */
    9,      /* colNo */
    "sqrt", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\elfun\\sqrt.m" /* pName
                                                                       */
};

/* Function Definitions */
real_T custom_power_iteration_noops(
    const emlrtStack *sp, const real_T Ahat_1[686], const real_T Ahat_2[686],
    const real_T Bhat_minus[196], const real_T Bhat_plus[196],
    const real_T Shat[98], real_T lx1_inv, real_T lx2_inv, real_T x[105],
    real_T xi[105], real_T u[30], real_T s)
{
  __m128d r;
  __m128d r1;
  __m128d r2;
  emlrtStack st;
  real_T Z[98];
  real_T b_Z[98];
  real_T varargin_2[98];
  real_T w_k[98];
  real_T b_Shat[91];
  real_T b_varargin_2[28];
  real_T a[15];
  real_T b_y[15];
  real_T c_y[15];
  real_T y[15];
  real_T b_s;
  real_T sigma;
  int32_T b_i;
  int32_T i;
  int32_T ii;
  int32_T j;
  int32_T k;
  int32_T s_tmp;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  emlrtMEXProfilingFunctionEntry((char_T *)d_custom_power_iteration_noops_,
                                 isMexOutdated);
  /* CUSTOM_POWER_ITERATION Computes max eigenvalue squared of dynamics matrix
   */
  /*    Input variables are problem variables obtained after L2-hypersphere */
  /*    preconditioning */
  emlrtMEXProfilingStatement(14, isMexOutdated);
  emlrtMEXProfilingStatement(15, isMexOutdated);
  emlrtMEXProfilingStatement(16, isMexOutdated);
  emlrtMEXProfilingStatement(17, isMexOutdated);
  emlrtMEXProfilingStatement(19, isMexOutdated);
  vecnorm(x, a);
  for (k = 0; k <= 12; k += 2) {
    r = _mm_loadu_pd(&a[k]);
    _mm_storeu_pd(&y[k], _mm_mul_pd(r, r));
  }
  y[14] = a[14] * a[14];
  vecnorm(xi, a);
  for (k = 0; k <= 12; k += 2) {
    r = _mm_loadu_pd(&a[k]);
    _mm_storeu_pd(&b_y[k], _mm_mul_pd(r, r));
  }
  b_y[14] = a[14] * a[14];
  b_vecnorm(u, a);
  for (k = 0; k <= 12; k += 2) {
    r = _mm_loadu_pd(&a[k]);
    r1 = _mm_loadu_pd(&y[k]);
    r2 = _mm_loadu_pd(&b_y[k]);
    _mm_storeu_pd(&c_y[k], _mm_add_pd(_mm_add_pd(r1, r2), _mm_mul_pd(r, r)));
  }
  c_y[14] = (y[14] + b_y[14]) + a[14] * a[14];
  sigma = s * s + sumColumnB(c_y);
  emlrtMEXProfilingStatement(20, isMexOutdated);
  st.site = &b_emlrtRSI;
  if (sigma < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  sigma = muDoubleScalarSqrt(sigma);
  emlrtMEXProfilingStatement(21, isMexOutdated);
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j < 5000)) {
    real_T d_y[637];
    real_T e_y[182];
    real_T f_y[182];
    real_T c_Z[26];
    real_T d_Z[13];
    real_T b_a;
    emlrtMEXProfilingStatement(22, isMexOutdated);
    emlrtMEXProfilingStatement(23, isMexOutdated);
    emlrtMEXProfilingStatement(24, isMexOutdated);
    emlrtMEXProfilingStatement(25, isMexOutdated);
    emlrtMEXProfilingStatement(26, isMexOutdated);
    emlrtMEXProfilingStatement(27, isMexOutdated);
    b_a = 1.0 / sigma;
    memcpy(&varargin_2[0], &x[0], 98U * sizeof(real_T));
    for (ii = 0; ii < 14; ii++) {
      for (i = 0; i < 7; i++) {
        b_s = 0.0;
        for (k = 0; k < 7; k++) {
          b_i = k * 7 + i;
          b_s += Ahat_1[(b_i % 7 + 7 * (b_i / 7)) + 49 * ii] *
                 varargin_2[k + 7 * ii];
        }
        w_k[i + 7 * ii] = b_s;
      }
    }
    memcpy(&varargin_2[0], &xi[0], 98U * sizeof(real_T));
    for (ii = 0; ii < 14; ii++) {
      for (i = 0; i < 7; i++) {
        b_s = 0.0;
        for (k = 0; k < 7; k++) {
          b_i = k * 7 + i;
          b_s += Ahat_2[(b_i % 7 + 7 * (b_i / 7)) + 49 * ii] *
                 varargin_2[k + 7 * ii];
        }
        Z[i + 7 * ii] = b_s;
      }
    }
    memcpy(&b_varargin_2[0], &u[0], 28U * sizeof(real_T));
    for (ii = 0; ii < 14; ii++) {
      s_tmp = ii << 1;
      for (i = 0; i < 7; i++) {
        varargin_2[i + 7 * ii] =
            Bhat_minus[i % 7 + 14 * ii] * b_varargin_2[s_tmp] +
            Bhat_minus[((i + 7) % 7 + 14 * ii) + 7] * b_varargin_2[s_tmp + 1];
      }
    }
    memcpy(&b_varargin_2[0], &u[2], 28U * sizeof(real_T));
    for (ii = 0; ii < 14; ii++) {
      s_tmp = ii << 1;
      for (i = 0; i < 7; i++) {
        b_Z[i + 7 * ii] =
            Bhat_plus[i % 7 + 14 * ii] * b_varargin_2[s_tmp] +
            Bhat_plus[((i + 7) % 7 + 14 * ii) + 7] * b_varargin_2[s_tmp + 1];
      }
    }
    for (b_i = 0; b_i <= 96; b_i += 2) {
      __m128d r3;
      __m128d r4;
      __m128d r5;
      r = _mm_loadu_pd(&w_k[b_i]);
      r1 = _mm_loadu_pd(&Z[b_i]);
      r2 = _mm_loadu_pd(&varargin_2[b_i]);
      r3 = _mm_loadu_pd(&b_Z[b_i]);
      r4 = _mm_loadu_pd(&x[b_i + 7]);
      r5 = _mm_loadu_pd(&xi[b_i + 7]);
      _mm_storeu_pd(
          &w_k[b_i],
          _mm_mul_pd(
              _mm_set1_pd(b_a),
              _mm_sub_pd(
                  _mm_sub_pd(
                      _mm_add_pd(
                          _mm_add_pd(_mm_add_pd(_mm_add_pd(r, r1), r2), r3),
                          _mm_mul_pd(_mm_loadu_pd(&Shat[b_i]), _mm_set1_pd(s))),
                      _mm_mul_pd(_mm_set1_pd(lx1_inv), r4)),
                  _mm_mul_pd(_mm_set1_pd(lx2_inv), r5))));
    }
    emlrtMEXProfilingStatement(28, isMexOutdated);
    for (b_i = 0; b_i < 7; b_i++) {
      b_s = 0.0;
      for (s_tmp = 0; s_tmp < 7; s_tmp++) {
        b_s += Ahat_1[s_tmp + 7 * b_i] * w_k[s_tmp];
      }
      x[b_i] = b_s;
    }
    emlrtMEXProfilingStatement(29, isMexOutdated);
    for (b_i = 0; b_i < 7; b_i++) {
      b_s = 0.0;
      for (s_tmp = 0; s_tmp < 7; s_tmp++) {
        b_s += Ahat_2[s_tmp + 7 * b_i] * w_k[s_tmp];
      }
      xi[b_i] = b_s;
    }
    emlrtMEXProfilingStatement(30, isMexOutdated);
    for (b_i = 0; b_i < 2; b_i++) {
      b_s = 0.0;
      for (s_tmp = 0; s_tmp < 7; s_tmp++) {
        b_s += Bhat_minus[s_tmp + 7 * b_i] * w_k[s_tmp];
      }
      u[b_i] = b_s;
    }
    emlrtMEXProfilingStatement(31, isMexOutdated);
    emlrtMEXProfilingStatement(32, isMexOutdated);
    for (k = 0; k < 13; k++) {
      for (ii = 0; ii < 7; ii++) {
        for (i = 0; i < 7; i++) {
          d_y[(ii + 7 * i) + 49 * k] = Ahat_1[(i + 7 * ii) + 49 * (k + 1)];
        }
      }
    }
    for (ii = 0; ii < 13; ii++) {
      for (i = 0; i < 7; i++) {
        s = 0.0;
        for (k = 0; k < 7; k++) {
          b_i = k * 7 + i;
          s += d_y[(b_i % 7 + 7 * (b_i / 7)) + 49 * ii] * w_k[k + 7 * (ii + 1)];
        }
        x[i + 7 * (ii + 1)] = s - lx1_inv * w_k[i + 7 * ii];
      }
    }
    emlrtMEXProfilingStatement(33, isMexOutdated);
    for (k = 0; k < 13; k++) {
      for (ii = 0; ii < 7; ii++) {
        for (i = 0; i < 7; i++) {
          d_y[(ii + 7 * i) + 49 * k] = Ahat_2[(i + 7 * ii) + 49 * (k + 1)];
        }
      }
    }
    for (ii = 0; ii < 13; ii++) {
      for (i = 0; i < 7; i++) {
        s = 0.0;
        for (k = 0; k < 7; k++) {
          b_i = k * 7 + i;
          s += d_y[(b_i % 7 + 7 * (b_i / 7)) + 49 * ii] * w_k[k + 7 * (ii + 1)];
        }
        xi[i + 7 * (ii + 1)] = s - lx2_inv * w_k[i + 7 * ii];
      }
    }
    emlrtMEXProfilingStatement(34, isMexOutdated);
    for (k = 0; k < 13; k++) {
      for (ii = 0; ii < 2; ii++) {
        for (i = 0; i < 7; i++) {
          e_y[(ii + (i << 1)) + 14 * k] =
              Bhat_minus[(i + 7 * ii) + 14 * (k + 1)];
        }
      }
    }
    for (k = 0; k < 13; k++) {
      for (ii = 0; ii < 2; ii++) {
        s = 0.0;
        for (i = 0; i < 7; i++) {
          s_tmp = ii + (i << 1);
          f_y[s_tmp + 14 * k] = Bhat_plus[(i + 7 * ii) + 14 * k];
          s += e_y[(s_tmp % 2 + ((s_tmp / 2) << 1)) + 14 * k] *
               w_k[i + 7 * (k + 1)];
        }
        c_Z[ii + (k << 1)] = s;
      }
    }
    for (ii = 0; ii < 13; ii++) {
      for (i = 0; i < 2; i++) {
        s = 0.0;
        for (k = 0; k < 7; k++) {
          b_i = (k << 1) + i;
          s += f_y[(b_i % 2 + ((b_i / 2) << 1)) + 14 * ii] * w_k[k + 7 * ii];
        }
        u[i + ((ii + 1) << 1)] = c_Z[i + (ii << 1)] + s;
      }
    }
    emlrtMEXProfilingStatement(35, isMexOutdated);
    memcpy(&b_Shat[0], &Shat[7], 91U * sizeof(real_T));
    for (ii = 0; ii < 13; ii++) {
      b_s = 0.0;
      for (k = 0; k < 7; k++) {
        b_s += b_Shat[k + 7 * ii] * w_k[k + 7 * (ii + 1)];
      }
      d_Z[ii] = b_s;
    }
    b_s = d_Z[0];
    for (k = 0; k < 12; k++) {
      b_s += d_Z[k + 1];
    }
    b_a = 0.0;
    for (b_i = 0; b_i < 7; b_i++) {
      b_a += Shat[b_i] * w_k[b_i];
    }
    s = b_a + b_s;
    emlrtMEXProfilingStatement(36, isMexOutdated);
    r = _mm_loadu_pd(&w_k[91]);
    r1 = _mm_set1_pd(-lx1_inv);
    _mm_storeu_pd(&x[98], _mm_mul_pd(r1, r));
    r = _mm_loadu_pd(&w_k[93]);
    _mm_storeu_pd(&x[100], _mm_mul_pd(r1, r));
    r = _mm_loadu_pd(&w_k[95]);
    _mm_storeu_pd(&x[102], _mm_mul_pd(r1, r));
    x[104] = -lx1_inv * w_k[97];
    emlrtMEXProfilingStatement(37, isMexOutdated);
    r = _mm_loadu_pd(&w_k[91]);
    r1 = _mm_set1_pd(-lx2_inv);
    _mm_storeu_pd(&xi[98], _mm_mul_pd(r1, r));
    r = _mm_loadu_pd(&w_k[93]);
    _mm_storeu_pd(&xi[100], _mm_mul_pd(r1, r));
    r = _mm_loadu_pd(&w_k[95]);
    _mm_storeu_pd(&xi[102], _mm_mul_pd(r1, r));
    xi[104] = -lx2_inv * w_k[97];
    emlrtMEXProfilingStatement(38, isMexOutdated);
    for (b_i = 0; b_i < 2; b_i++) {
      b_s = 0.0;
      for (s_tmp = 0; s_tmp < 7; s_tmp++) {
        b_s += Bhat_plus[(s_tmp + 7 * b_i) + 182] * w_k[s_tmp + 91];
      }
      u[b_i + 28] = b_s;
    }
    emlrtMEXProfilingStatement(39, isMexOutdated);
    vecnorm(x, a);
    for (k = 0; k <= 12; k += 2) {
      r = _mm_loadu_pd(&a[k]);
      _mm_storeu_pd(&y[k], _mm_mul_pd(r, r));
    }
    y[14] = a[14] * a[14];
    vecnorm(xi, a);
    for (k = 0; k <= 12; k += 2) {
      r = _mm_loadu_pd(&a[k]);
      _mm_storeu_pd(&b_y[k], _mm_mul_pd(r, r));
    }
    b_y[14] = a[14] * a[14];
    b_vecnorm(u, a);
    for (k = 0; k <= 12; k += 2) {
      r = _mm_loadu_pd(&a[k]);
      r1 = _mm_loadu_pd(&y[k]);
      r2 = _mm_loadu_pd(&b_y[k]);
      _mm_storeu_pd(&c_y[k], _mm_add_pd(_mm_add_pd(r1, r2), _mm_mul_pd(r, r)));
    }
    c_y[14] = (y[14] + b_y[14]) + a[14] * a[14];
    b_s = s * s + sumColumnB(c_y);
    emlrtMEXProfilingStatement(40, isMexOutdated);
    st.site = &d_emlrtRSI;
    if (b_s < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    b_s = muDoubleScalarSqrt(b_s);
    emlrtMEXProfilingStatement(41, isMexOutdated);
    if (muDoubleScalarAbs(b_s - sigma) <=
        0.001 * muDoubleScalarMax(b_s, sigma) + 0.001) {
      emlrtMEXProfilingStatement(42, isMexOutdated);
      exitg1 = true;
    } else {
      emlrtMEXProfilingStatement(43, isMexOutdated);
      if (j + 1 < 5000) {
        emlrtMEXProfilingStatement(44, isMexOutdated);
        sigma = b_s;
        emlrtMEXProfilingStatement(45, isMexOutdated);
      }
      emlrtMEXProfilingStatement(46, isMexOutdated);
      j++;
    }
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  emlrtMEXProfilingStatement(47, isMexOutdated);
  sigma = 1.1 * b_s;
  emlrtMEXProfilingStatement(48, isMexOutdated);
  emlrtMEXProfilingFunctionExit(isMexOutdated);
  return sigma;
}

/* End of code generation (custom_power_iteration_noops.c) */
