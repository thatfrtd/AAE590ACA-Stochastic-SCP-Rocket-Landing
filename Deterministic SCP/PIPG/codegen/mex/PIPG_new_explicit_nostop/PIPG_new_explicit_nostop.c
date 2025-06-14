/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_new_explicit_nostop.c
 *
 * Code generation for function 'PIPG_new_explicit_nostop'
 *
 */

/* Include files */
#include "PIPG_new_explicit_nostop.h"
#include "PIPG_new_explicit_nostop_data.h"
#include "PIPG_new_explicit_nostop_emxutil.h"
#include "PIPG_new_explicit_nostop_types.h"
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "tic.h"
#include "toc.h"
#include "blas.h"
#include "emlrt.h"
#include "mwmathutil.h"
#include <emmintrin.h>
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    20,                         /* lineNo */
    "PIPG_new_explicit_nostop", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    33,                         /* lineNo */
    "PIPG_new_explicit_nostop", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    77,                         /* lineNo */
    "PIPG_new_explicit_nostop", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m" /* pathName */
};

static emlrtDCInfo emlrtDCI = {
    16,                         /* lineNo */
    36,                         /* colNo */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m", /* pName */
    1                                  /* checkKind */
};

static emlrtBCInfo emlrtBCI = {
    1,                          /* iFirst */
    1,                          /* iLast */
    16,                         /* lineNo */
    36,                         /* colNo */
    "[size(z_ref, j_max)]",     /* aName */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m", /* pName */
    0                                  /* checkKind */
};

static emlrtRTEInfo emlrtRTEI = {
    44,                         /* lineNo */
    9,                          /* colNo */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m" /* pName */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,                         /* iFirst */
    -1,                         /* iLast */
    50,                         /* lineNo */
    13,                         /* colNo */
    "zhis",                     /* aName */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m", /* pName */
    0                                  /* checkKind */
};

static emlrtRTEInfo b_emlrtRTEI = {
    13,     /* lineNo */
    9,      /* colNo */
    "sqrt", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\elfun\\sqrt.m" /* pName
                                                                       */
};

static emlrtDCInfo b_emlrtDCI = {
    14,                         /* lineNo */
    30,                         /* colNo */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m", /* pName */
    4                                  /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = {
    14,                         /* lineNo */
    30,                         /* colNo */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m", /* pName */
    1                                  /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = {
    15,                         /* lineNo */
    30,                         /* colNo */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m", /* pName */
    1                                  /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = {
    42,                         /* lineNo */
    14,                         /* colNo */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m", /* pName */
    1                                  /* checkKind */
};

static emlrtRTEInfo d_emlrtRTEI = {
    14,                         /* lineNo */
    1,                          /* colNo */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m" /* pName */
};

static emlrtRTEInfo e_emlrtRTEI = {
    15,                         /* lineNo */
    1,                          /* colNo */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m" /* pName */
};

static emlrtRTEInfo f_emlrtRTEI = {
    42,                         /* lineNo */
    1,                          /* colNo */
    "PIPG_new_explicit_nostop", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit_nostop.m" /* pName */
};

/* Function Definitions */
void PIPG_new_explicit_nostop(
    const emlrtStack *sp, const real_T qhat[1000], const real_T Hhat[3000],
    const real_T h[3], const real_T L[1000000], const real_T L_inv[1000000],
    real_T lambda, real_T sigma, real_T omega, real_T rho, real_T tol_abs,
    real_T tol_rel, real_T tol_infeas, real_T j_check, real_T j_max,
    const real_T z_ref[1000], const real_T what_ref[3], real_T z_star[1000],
    real_T what_star[3], struct0_T *sol_info)
{
  static const char_T t0_Value[8] = {'U', 'n', 's', 'o', 'l', 'v', 'e', 'd'};
  ptrdiff_t b_k_t;
  ptrdiff_t b_lda_t;
  ptrdiff_t b_ldb_t;
  ptrdiff_t b_ldc_t;
  ptrdiff_t b_m_t;
  ptrdiff_t b_n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  emlrtStack st;
  emlrtTimespec expl_temp;
  real_T b_zeta_jp1[1000];
  real_T y[1000];
  real_T zeta_jp1[1000];
  real_T zhat_j[1000];
  real_T zhat_jp1[1000];
  real_T eta_jp1[3];
  real_T alpha;
  real_T alpha1 = 0.0;
  real_T b_alpha1 = 0.0;
  real_T b_beta1 = 0.0;
  real_T beta;
  real_T beta1 = 0.0;
  real_T d;
  real_T d1;
  int32_T b_i;
  int32_T i;
  int32_T j;
  int32_T loop_ub_tmp;
  char_T TRANSA1 = '\x00';
  char_T TRANSB1 = '\x00';
  char_T b_TRANSA1 = '\x00';
  char_T b_TRANSB1 = '\x00';
  (void)tol_abs;
  (void)tol_rel;
  (void)tol_infeas;
  (void)j_check;
  st.prev = sp;
  st.tls = sp->tls;
  emlrtMEXProfilingFunctionEntry((char_T *)c_PIPG_new_explicit_nostop_comp,
                                 isMexOutdated);
  /* PIPG Proportional Integral Projected Gradient convex conic solver */
  /*    Extrapolated Proportional Integral Projected Gradient (xPIPG) */
  /*  rho in [1.5, 1.9] is usually a good choice - should converge about twice
   */
  /*  as fast as not doing extrapolation (rho = 1) */
  /*  sigma should be >= norm(H) */
  /*  lambda should be >= norm(Q) */
  emlrtMEXProfilingStatement(2, isMexOutdated);
  sol_info->primal_infeasible = false;
  emlrtMEXProfilingStatement(3, isMexOutdated);
  sol_info->dual_infeasible = false;
  emlrtMEXProfilingStatement(4, isMexOutdated);
  for (i = 0; i < 8; i++) {
    sol_info->solution_status.Value[i] = t0_Value[i];
  }
  emlrtMEXProfilingStatement(5, isMexOutdated);
  sol_info->iterations = -1.0;
  emlrtMEXProfilingStatement(6, isMexOutdated);
  if (!(j_max >= 0.0)) {
    emlrtNonNegativeCheckR2012b(j_max, &b_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = (int32_T)muDoubleScalarFloor(j_max);
  if (j_max != i) {
    emlrtIntegerCheckR2012b(j_max, &c_emlrtDCI, (emlrtConstCTX)sp);
  }
  loop_ub_tmp = (int32_T)j_max;
  b_i = sol_info->zhat_inf_dj->size[0];
  sol_info->zhat_inf_dj->size[0] = loop_ub_tmp;
  emxEnsureCapacity_real_T(sp, sol_info->zhat_inf_dj, b_i, &d_emlrtRTEI);
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    sol_info->zhat_inf_dj->data[b_i] = 0.0;
  }
  emlrtMEXProfilingStatement(7, isMexOutdated);
  if (loop_ub_tmp != i) {
    emlrtIntegerCheckR2012b(j_max, &d_emlrtDCI, (emlrtConstCTX)sp);
  }
  b_i = sol_info->what_inf_dj->size[0];
  sol_info->what_inf_dj->size[0] = loop_ub_tmp;
  emxEnsureCapacity_real_T(sp, sol_info->what_inf_dj, b_i, &e_emlrtRTEI);
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    sol_info->what_inf_dj->data[b_i] = 0.0;
  }
  emlrtMEXProfilingStatement(8, isMexOutdated);
  if (loop_ub_tmp != i) {
    emlrtIntegerCheckR2012b(j_max, &emlrtDCI, (emlrtConstCTX)sp);
  }
  if (j_max < 1.0) {
    emlrtDynamicBoundsCheckR2012b(0, 1, 1, &emlrtBCI, (emlrtConstCTX)sp);
  }
  emlrtMEXProfilingStatement(10, isMexOutdated);
  st.site = &emlrtRSI;
  expl_temp = tic(&st);
  /*  Transform previous primal solution */
  emlrtMEXProfilingStatement(11, isMexOutdated);
  mtimes(L, z_ref, zeta_jp1);
  /*  Initialize transformed primal variables */
  emlrtMEXProfilingStatement(12, isMexOutdated);
  /*  Intialize transformed dual variable */
  emlrtMEXProfilingStatement(13, isMexOutdated);
  eta_jp1[0] = what_ref[0];
  eta_jp1[1] = what_ref[1];
  eta_jp1[2] = what_ref[2];
  /*  Calculate step-size (sigma from power iteration) */
  emlrtMEXProfilingStatement(15, isMexOutdated);
  st.site = &c_emlrtRSI;
  alpha = lambda * lambda + 4.0 * omega * sigma;
  if (alpha < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &b_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  alpha = muDoubleScalarSqrt(alpha);
  alpha = 2.0 / (lambda + alpha);
  emlrtMEXProfilingStatement(16, isMexOutdated);
  beta = omega * alpha;
  /*  Initialize */
  emlrtMEXProfilingStatement(17, isMexOutdated);
  memcpy(&zhat_j[0], &zeta_jp1[0], 1000U * sizeof(real_T));
  emlrtMEXProfilingStatement(18, isMexOutdated);
  memset(&zhat_jp1[0], 0, 1000U * sizeof(real_T));
  emlrtMEXProfilingStatement(19, isMexOutdated);
  what_star[0] = 0.0;
  what_star[1] = 0.0;
  what_star[2] = 0.0;
  emlrtMEXProfilingStatement(20, isMexOutdated);
  if ((int32_T)j_max != i) {
    emlrtIntegerCheckR2012b(j_max, &e_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = sol_info->zhis->size[0] * sol_info->zhis->size[1] *
      sol_info->zhis->size[2];
  sol_info->zhis->size[0] = 1000;
  sol_info->zhis->size[1] = 1;
  sol_info->zhis->size[2] = loop_ub_tmp;
  emxEnsureCapacity_real_T(sp, sol_info->zhis, i, &f_emlrtRTEI);
  b_i = 1000 * (int32_T)j_max;
  for (i = 0; i < b_i; i++) {
    sol_info->zhis->data[i] = 0.0;
  }
  emlrtMEXProfilingStatement(21, isMexOutdated);
  emlrtForLoopVectorCheckR2021a(1.0, 1.0, j_max, mxDOUBLE_CLASS, (int32_T)j_max,
                                &emlrtRTEI, (emlrtConstCTX)sp);
  if (j_max - 1.0 >= 0.0) {
    TRANSB1 = 'N';
    TRANSA1 = 'T';
    alpha1 = 1.0;
    m_t = (ptrdiff_t)1000;
    n_t = (ptrdiff_t)1;
    k_t = (ptrdiff_t)3;
    lda_t = (ptrdiff_t)3;
    ldb_t = (ptrdiff_t)3;
    ldc_t = (ptrdiff_t)1000;
    b_TRANSB1 = 'N';
    b_TRANSA1 = 'N';
    b_alpha1 = 1.0;
    b_m_t = (ptrdiff_t)3;
    b_n_t = (ptrdiff_t)1;
    b_k_t = (ptrdiff_t)1000;
    b_lda_t = (ptrdiff_t)3;
    b_ldb_t = (ptrdiff_t)1000;
    b_ldc_t = (ptrdiff_t)3;
    d = 1.0 - rho;
    d1 = 1.0 - rho;
  }
  if (loop_ub_tmp - 1 >= 0) {
  }
  for (j = 0; j < loop_ub_tmp; j++) {
    __m128d r;
    __m128d r1;
    real_T absxk;
    real_T scale;
    real_T t;
    real_T varargin_1;
    /*     %% Projected gradient step */
    emlrtMEXProfilingStatement(22, isMexOutdated);
    dgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &alpha1, (real_T *)&Hhat[0],
          &lda_t, &eta_jp1[0], &ldb_t, &beta1, &y[0], &ldc_t);
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zeta_jp1[i]);
      r1 = _mm_loadu_pd(&y[i]);
      _mm_storeu_pd(
          &b_zeta_jp1[i],
          _mm_sub_pd(
              r, _mm_mul_pd(
                     _mm_set1_pd(alpha),
                     _mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_set1_pd(lambda), r),
                                           _mm_loadu_pd(&qhat[i])),
                                r1))));
    }
    mtimes(L_inv, b_zeta_jp1, zhat_jp1);
    emlrtMEXProfilingFunctionEntry((char_T *)ball_project_complete_name,
                                   isMexOutdated);
    /* PROJECT Summary of this method goes here */
    /*    Detailed explanation goes here */
    emlrtMEXProfilingStatement(1, isMexOutdated);
    emlrtMEXProfilingStatement(3, isMexOutdated);
    emlrtMEXProfilingStatement(4, isMexOutdated);
    scale = 3.3121686421112381E-170;
    absxk = muDoubleScalarAbs(zhat_jp1[0]);
    if (absxk > 3.3121686421112381E-170) {
      varargin_1 = 1.0;
      scale = absxk;
    } else {
      t = absxk / 3.3121686421112381E-170;
      varargin_1 = t * t;
    }
    absxk = muDoubleScalarAbs(zhat_jp1[1]);
    if (absxk > scale) {
      t = scale / absxk;
      varargin_1 = varargin_1 * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      varargin_1 += t * t;
    }
    r = _mm_loadu_pd(&zhat_jp1[0]);
    _mm_storeu_pd(
        &zhat_jp1[0],
        _mm_mul_pd(
            _mm_set1_pd(
                16.0 / muDoubleScalarMax(scale * muDoubleScalarSqrt(varargin_1),
                                         16.0)),
            r));
    emlrtMEXProfilingStatement(5, isMexOutdated);
    emlrtMEXProfilingStatement(6, isMexOutdated);
    emlrtMEXProfilingFunctionExit(isMexOutdated);
    emlrtMEXProfilingStatement(23, isMexOutdated);
    emlrtMEXProfilingFunctionEntry((char_T *)singleton_project_complete_name,
                                   isMexOutdated);
    /* PROJECT Summary of this method goes here */
    /*    Detailed explanation goes here */
    emlrtMEXProfilingStatement(1, isMexOutdated);
    memcpy(&y[0], &zhat_jp1[0], 1000U * sizeof(real_T));
    emlrtMEXProfilingStatement(2, isMexOutdated);
    y[9] = 1.0;
    y[10] = 2.0;
    y[11] = 3.0;
    emlrtMEXProfilingStatement(3, isMexOutdated);
    emlrtMEXProfilingFunctionExit(isMexOutdated);
    mtimes(L, y, zhat_jp1);
    /* aff_violation(j) = norm(Hhat * zhat_jp1 - h); */
    emlrtMEXProfilingStatement(24, isMexOutdated);
    if ((int32_T)((uint32_T)j + 1U) > sol_info->zhis->size[2]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)j + 1U), 1,
                                    sol_info->zhis->size[2], &b_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    for (i = 0; i < 1000; i++) {
      sol_info->zhis->data[i + 1000 * j] = zhat_j[i];
    }
    /*     %% Proportional-integral feedback of affine equality constraint
     * violation */
    emlrtMEXProfilingStatement(25, isMexOutdated);
    for (b_i = 0; b_i <= 998; b_i += 2) {
      r = _mm_loadu_pd(&zhat_jp1[b_i]);
      r1 = _mm_loadu_pd(&zeta_jp1[b_i]);
      _mm_storeu_pd(&y[b_i], _mm_sub_pd(_mm_mul_pd(_mm_set1_pd(2.0), r), r1));
    }
    dgemm(&b_TRANSA1, &b_TRANSB1, &b_m_t, &b_n_t, &b_k_t, &b_alpha1,
          (real_T *)&Hhat[0], &b_lda_t, &y[0], &b_ldb_t, &b_beta1,
          &what_star[0], &b_ldc_t);
    r = _mm_loadu_pd(&what_star[0]);
    r1 = _mm_loadu_pd(&eta_jp1[0]);
    _mm_storeu_pd(
        &what_star[0],
        _mm_add_pd(r1, _mm_mul_pd(_mm_set1_pd(beta),
                                  _mm_sub_pd(r, _mm_loadu_pd(&h[0])))));
    what_star[2] = eta_jp1[2] + beta * (what_star[2] - h[2]);
    /*     %% Extrapolate transformed primal variables */
    emlrtMEXProfilingStatement(26, isMexOutdated);
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zeta_jp1[i]);
      r1 = _mm_loadu_pd(&zhat_jp1[i]);
      _mm_storeu_pd(&zeta_jp1[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd(d), r),
                                             _mm_mul_pd(_mm_set1_pd(rho), r1)));
    }
    /*     %% Extrapolate transformed dual variable */
    emlrtMEXProfilingStatement(27, isMexOutdated);
    r = _mm_loadu_pd(&eta_jp1[0]);
    r1 = _mm_loadu_pd(&what_star[0]);
    _mm_storeu_pd(&eta_jp1[0], _mm_add_pd(_mm_mul_pd(_mm_set1_pd(d1), r),
                                          _mm_mul_pd(_mm_set1_pd(rho), r1)));
    eta_jp1[2] = d1 * eta_jp1[2] + rho * what_star[2];
    /*     %% Update values for next iteration */
    emlrtMEXProfilingStatement(28, isMexOutdated);
    memcpy(&zhat_j[0], &zhat_jp1[0], 1000U * sizeof(real_T));
    emlrtMEXProfilingStatement(30, isMexOutdated);
    emlrtMEXProfilingStatement(31, isMexOutdated);
    emlrtMEXProfilingStatement(32, isMexOutdated);
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  /*  Recover original primal variables */
  emlrtMEXProfilingStatement(33, isMexOutdated);
  mtimes(L_inv, zhat_jp1, z_star);
  /*  Retain transformed dual variable */
  emlrtMEXProfilingStatement(34, isMexOutdated);
  emlrtMEXProfilingStatement(35, isMexOutdated);
  st.site = &h_emlrtRSI;
  alpha = toc(&st, expl_temp.tv_sec, expl_temp.tv_nsec);
  emlrtMEXProfilingStatement(36, isMexOutdated);
  sol_info->time = alpha;
  /*  %  */
  /* fprintf("PIPG Time: %.3f ms\n", t2 * 1000 ) */
  /*   */
  /*  figure */
  /*  plot(aff_violation) */
  /*  yscale("log") */
  /*   */
  /*  figure */
  /*  zhis = L_inv * zhis; */
  /*  scatter(zhis(25 * 7 + 1, :), zhis(25 * 7 + 2, :)) */
  /*  axis equal */
  /*  grid on */
  /*   */
  /*  figure */
  /*  zhis = L_inv * zhis; */
  /*  plot(vecnorm(zhis(1:2, :) - zhis(1:2, end))) */
  /*  yscale("log") */
  /*  grid on */
  /*  %  */
  emlrtMEXProfilingStatement(37, isMexOutdated);
  emlrtMEXProfilingStatement(38, isMexOutdated);
  emlrtMEXProfilingFunctionExit(isMexOutdated);
}

/* End of code generation (PIPG_new_explicit_nostop.c) */
