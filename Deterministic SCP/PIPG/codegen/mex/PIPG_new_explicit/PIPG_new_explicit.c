/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_new_explicit.c
 *
 * Code generation for function 'PIPG_new_explicit'
 *
 */

/* Include files */
#include "PIPG_new_explicit.h"
#include "PIPG_new_explicit_data.h"
#include "PIPG_new_explicit_emxutil.h"
#include "PIPG_new_explicit_mexutil.h"
#include "PIPG_new_explicit_types.h"
#include "mod.h"
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "stopping.h"
#include "tic.h"
#include "toc.h"
#include "vecnorm.h"
#include "emlrt.h"
#include "mwmathutil.h"
#include <emmintrin.h>
#include <stdio.h>
#include <string.h>

/* Variable Definitions */
static emlrtRTEInfo emlrtRTEI = {
    114,                 /* lineNo */
    1,                   /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pName */
};

static emlrtRSInfo emlrtRSI = {
    21,                  /* lineNo */
    "PIPG_new_explicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    24,                  /* lineNo */
    "PIPG_new_explicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    34,                  /* lineNo */
    "PIPG_new_explicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    47,                  /* lineNo */
    "PIPG_new_explicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    48,                  /* lineNo */
    "PIPG_new_explicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    54,                  /* lineNo */
    "PIPG_new_explicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    106,                 /* lineNo */
    "PIPG_new_explicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    111,                 /* lineNo */
    "PIPG_new_explicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI = {
    114,                 /* lineNo */
    "PIPG_new_explicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pathName */
};

static emlrtRSInfo n_emlrtRSI =
    {
        94,                  /* lineNo */
        "eml_mtimes_helper", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_"
        "helper.m" /* pathName */
};

static emlrtRSInfo q_emlrtRSI = {
    44,       /* lineNo */
    "mpower", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\matfun\\mpower.m" /* pathName
                                                                          */
};

static emlrtRSInfo r_emlrtRSI =
    {
        71,      /* lineNo */
        "power", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\power.m" /* pathName
                                                                          */
};

static emlrtRSInfo t_emlrtRSI = {
    38,        /* lineNo */
    "fprintf", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\iofun\\fprintf.m" /* pathName
                                                                          */
};

static emlrtMCInfo emlrtMCI = {
    23,              /* lineNo */
    13,              /* colNo */
    "Ballz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Ballz.m" /* pName */
};

static emlrtMCInfo b_emlrtMCI = {
    66,        /* lineNo */
    18,        /* colNo */
    "fprintf", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\iofun\\fprintf.m" /* pName
                                                                          */
};

static emlrtDCInfo emlrtDCI = {
    17,                  /* lineNo */
    36,                  /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    1                           /* checkKind */
};

static emlrtBCInfo emlrtBCI = {
    1,                      /* iFirst */
    1,                      /* iLast */
    17,                     /* lineNo */
    36,                     /* colNo */
    "[size(z_ref, j_max)]", /* aName */
    "PIPG_new_explicit",    /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    0                           /* checkKind */
};

static emlrtRTEInfo b_emlrtRTEI = {
    45,                  /* lineNo */
    9,                   /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pName */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,                  /* iFirst */
    -1,                  /* iLast */
    51,                  /* lineNo */
    13,                  /* colNo */
    "zhis",              /* aName */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    0                           /* checkKind */
};

static emlrtRTEInfo c_emlrtRTEI = {
    13,     /* lineNo */
    9,      /* colNo */
    "sqrt", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\elfun\\sqrt.m" /* pName
                                                                       */
};

static emlrtDCInfo b_emlrtDCI = {
    15,                  /* lineNo */
    30,                  /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    4                           /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = {
    15,                  /* lineNo */
    30,                  /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    1                           /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = {
    16,                  /* lineNo */
    30,                  /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    1                           /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = {
    43,                  /* lineNo */
    14,                  /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    1                           /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,                     /* iFirst */
    -1,                     /* iLast */
    65,                     /* lineNo */
    30,                     /* colNo */
    "sol_info.zhat_inf_dj", /* aName */
    "PIPG_new_explicit",    /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    0                           /* checkKind */
};

static emlrtDCInfo f_emlrtDCI = {
    65,                  /* lineNo */
    30,                  /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    1                           /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,                     /* iFirst */
    -1,                     /* iLast */
    66,                     /* lineNo */
    30,                     /* colNo */
    "sol_info.what_inf_dj", /* aName */
    "PIPG_new_explicit",    /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    0                           /* checkKind */
};

static emlrtDCInfo g_emlrtDCI = {
    66,                  /* lineNo */
    30,                  /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m", /* pName */
    1                           /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = {
    1,                    /* iFirst */
    1000,                 /* iLast */
    35,                   /* lineNo */
    13,                   /* colNo */
    "x",                  /* aName */
    "Singletonz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Singletonz.m", /* pName */
    3                                  /* checkKind */
};

static emlrtDCInfo h_emlrtDCI = {
    35,                   /* lineNo */
    13,                   /* colNo */
    "Singletonz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Singletonz.m", /* pName */
    1                                  /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = {
    1,                    /* iFirst */
    1000,                 /* iLast */
    31,                   /* lineNo */
    28,                   /* colNo */
    "z",                  /* aName */
    "Singletonz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Singletonz.m", /* pName */
    0                                  /* checkKind */
};

static emlrtDCInfo i_emlrtDCI = {
    31,                   /* lineNo */
    28,                   /* colNo */
    "Singletonz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Singletonz.m", /* pName */
    1                                  /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = {
    1,                    /* iFirst */
    1000,                 /* iLast */
    30,                   /* lineNo */
    28,                   /* colNo */
    "z",                  /* aName */
    "Singletonz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Singletonz.m", /* pName */
    0                                  /* checkKind */
};

static emlrtDCInfo j_emlrtDCI = {
    30,                   /* lineNo */
    28,                   /* colNo */
    "Singletonz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Singletonz.m", /* pName */
    1                                  /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = {
    1,               /* iFirst */
    1000,            /* iLast */
    32,              /* lineNo */
    28,              /* colNo */
    "z",             /* aName */
    "Ballz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Ballz.m", /* pName */
    0                             /* checkKind */
};

static emlrtDCInfo k_emlrtDCI = {
    32,              /* lineNo */
    28,              /* colNo */
    "Ballz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Ballz.m", /* pName */
    1                             /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = {
    1,               /* iFirst */
    1000,            /* iLast */
    31,              /* lineNo */
    28,              /* colNo */
    "z",             /* aName */
    "Ballz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Ballz.m", /* pName */
    0                             /* checkKind */
};

static emlrtDCInfo l_emlrtDCI = {
    31,              /* lineNo */
    28,              /* colNo */
    "Ballz/project", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Ballz.m", /* pName */
    1                             /* checkKind */
};

static emlrtRTEInfo e_emlrtRTEI = {
    15,                  /* lineNo */
    1,                   /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pName */
};

static emlrtRTEInfo f_emlrtRTEI = {
    16,                  /* lineNo */
    1,                   /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pName */
};

static emlrtRTEInfo g_emlrtRTEI = {
    43,                  /* lineNo */
    1,                   /* colNo */
    "PIPG_new_explicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_new_explicit.m" /* pName */
};

static emlrtRSInfo u_emlrtRSI = {
    23,              /* lineNo */
    "Ballz/project", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Ballz.m" /* pathName */
};

static emlrtRSInfo v_emlrtRSI = {
    66,        /* lineNo */
    "fprintf", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\iofun\\fprintf.m" /* pathName
                                                                          */
};

/* Function Declarations */
static void b_error(const emlrtStack *sp, const mxArray *m,
                    emlrtMCInfo *location);

static const mxArray *e_emlrt_marshallOut(const emlrtStack *sp);

static const mxArray *feval(const emlrtStack *sp, const mxArray *m1,
                            const mxArray *m2, const mxArray *m3,
                            const mxArray *m4, emlrtMCInfo *location);

/* Function Definitions */
static void b_error(const emlrtStack *sp, const mxArray *m,
                    emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = m;
  emlrtCallMATLABR2012b((emlrtConstCTX)sp, 0, NULL, 1, &pArray, "error", true,
                        location);
}

static const mxArray *e_emlrt_marshallOut(const emlrtStack *sp)
{
  static const int32_T iv[2] = {1, 25};
  static const char_T u[25] = {'R', 'a', 'd', 'i', 'u', 's', ' ', 's', 'h',
                               'o', 'u', 'l', 'd', ' ', 'b', 'e', ' ', 'p',
                               'o', 's', 'i', 't', 'i', 'v', 'e'};
  const mxArray *b_y;
  const mxArray *m;
  const mxArray *propValues;
  const mxArray *y;
  const char_T *propClasses = "coder.internal.string";
  const char_T *propNames = "Value";
  y = NULL;
  m = NULL;
  emlrtAssign(
      &y, emlrtCreateClassInstance2022a((emlrtCTX)sp, "coder.internal.string"));
  m = NULL;
  b_y = NULL;
  propValues = emlrtCreateCharArray(2, &iv[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 25, propValues, &u[0]);
  emlrtAssign(&b_y, propValues);
  emlrtAssign(&m, b_y);
  propValues = m;
  emlrtSetAllProperties((emlrtCTX)sp, &y, 0, 1, (const char_T **)&propNames,
                        (const char_T **)&propClasses, &propValues);
  emlrtAssign(&y, emlrtConvertInstanceToRedirectSource(
                      (emlrtCTX)sp, y, 0, "coder.internal.string"));
  return y;
}

static const mxArray *feval(const emlrtStack *sp, const mxArray *m1,
                            const mxArray *m2, const mxArray *m3,
                            const mxArray *m4, emlrtMCInfo *location)
{
  const mxArray *pArrays[4];
  const mxArray *m;
  pArrays[0] = m1;
  pArrays[1] = m2;
  pArrays[2] = m3;
  pArrays[3] = m4;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m, 4, &pArrays[0],
                               "feval", true, location);
}

void PIPG_new_explicit(const emlrtStack *sp, const real_T qhat[1000],
                       const real_T Hhat[3000], const real_T h[3],
                       const Ballz *Dball, const Singletonz *Dsingleton,
                       const real_T L[1000000], const real_T L_inv[1000000],
                       real_T lambda, real_T sigma, real_T omega, real_T rho,
                       real_T tol_abs, real_T tol_rel, real_T tol_infeas,
                       real_T j_check, real_T j_max, const real_T z_ref[1000],
                       const real_T what_ref[3], real_T z_star[1000],
                       real_T what_star[3], struct0_T *sol_info)
{
  static const int32_T iv[2] = {1, 7};
  static const int32_T iv1[2] = {1, 20};
  static const char_T b_u[20] = {'P', 'I', 'P', 'G', ' ',  'T', 'i',
                                 'm', 'e', ':', ' ', '%',  '.', '3',
                                 'f', ' ', 'm', 's', '\\', 'n'};
  static const char_T cv2[10] = {'I', 'n', 'f', 'e', 'a',
                                 's', 'i', 'b', 'l', 'e'};
  static const char_T b_cv[8] = {'U', 'n', 's', 'o', 'l', 'v', 'e', 'd'};
  static const char_T cv1[7] = {'O', 'p', 't', 'i', 'm', 'a', 'l'};
  static const char_T u[7] = {'f', 'p', 'r', 'i', 'n', 't', 'f'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  emlrtTimespec expl_temp;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *m;
  const mxArray *y;
  real_T b_zhat_jp1[1000];
  real_T zeta_jp1[1000];
  real_T zhat_j[1000];
  real_T zhat_jp1[1000];
  real_T eta_jp1[3];
  real_T what_j[3];
  real_T alpha;
  real_T beta;
  real_T t2;
  real_T what_inf_dj;
  int32_T i;
  int32_T i1;
  int32_T loop_ub_tmp;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtMEXProfilingFunctionEntry((char_T *)PIPG_new_explicit_complete_name,
                                 isMexOutdated);
  /* PIPG Proportional Integral Projected Gradient convex conic solver */
  /*    Extrapolated Proportional Integral Projected Gradient (xPIPG) */
  /*  rho in [1.5, 1.9] is usually a good choice - should converge about twice
   */
  /*  as fast as not doing extrapolation (rho = 1) */
  /*  sigma should be >= norm(H) */
  /*  lambda should be >= norm(Q) */
  emlrtMEXProfilingStatement(3, isMexOutdated);
  sol_info->primal_infeasible = false;
  emlrtMEXProfilingStatement(4, isMexOutdated);
  sol_info->dual_infeasible = false;
  emlrtMEXProfilingStatement(5, isMexOutdated);
  sol_info->solution_status.Value.size[0] = 1;
  sol_info->solution_status.Value.size[1] = 8;
  for (i = 0; i < 8; i++) {
    sol_info->solution_status.Value.data[i] = b_cv[i];
  }
  emlrtMEXProfilingStatement(6, isMexOutdated);
  sol_info->iterations = -1.0;
  emlrtMEXProfilingStatement(7, isMexOutdated);
  if (!(j_max >= 0.0)) {
    emlrtNonNegativeCheckR2012b(j_max, &b_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = (int32_T)muDoubleScalarFloor(j_max);
  if (j_max != i) {
    emlrtIntegerCheckR2012b(j_max, &c_emlrtDCI, (emlrtConstCTX)sp);
  }
  loop_ub_tmp = (int32_T)j_max;
  i1 = sol_info->zhat_inf_dj->size[0];
  sol_info->zhat_inf_dj->size[0] = loop_ub_tmp;
  emxEnsureCapacity_real_T(sp, sol_info->zhat_inf_dj, i1, &e_emlrtRTEI);
  for (i1 = 0; i1 < loop_ub_tmp; i1++) {
    sol_info->zhat_inf_dj->data[i1] = 0.0;
  }
  emlrtMEXProfilingStatement(8, isMexOutdated);
  if (loop_ub_tmp != i) {
    emlrtIntegerCheckR2012b(j_max, &d_emlrtDCI, (emlrtConstCTX)sp);
  }
  i1 = sol_info->what_inf_dj->size[0];
  sol_info->what_inf_dj->size[0] = loop_ub_tmp;
  emxEnsureCapacity_real_T(sp, sol_info->what_inf_dj, i1, &f_emlrtRTEI);
  for (i1 = 0; i1 < loop_ub_tmp; i1++) {
    sol_info->what_inf_dj->data[i1] = 0.0;
  }
  emlrtMEXProfilingStatement(9, isMexOutdated);
  if (loop_ub_tmp != i) {
    emlrtIntegerCheckR2012b(j_max, &emlrtDCI, (emlrtConstCTX)sp);
  }
  if (j_max < 1.0) {
    emlrtDynamicBoundsCheckR2012b(0, 1, 1, &emlrtBCI, (emlrtConstCTX)sp);
  }
  emlrtMEXProfilingStatement(11, isMexOutdated);
  st.site = &emlrtRSI;
  expl_temp = tic(&st);
  /*  Transform previous primal solution */
  emlrtMEXProfilingStatement(12, isMexOutdated);
  st.site = &b_emlrtRSI;
  b_st.site = &n_emlrtRSI;
  mtimes(L, z_ref, zeta_jp1);
  /*  Initialize transformed primal variables */
  emlrtMEXProfilingStatement(13, isMexOutdated);
  /*  Intialize transformed dual variable */
  emlrtMEXProfilingStatement(14, isMexOutdated);
  eta_jp1[0] = what_ref[0];
  eta_jp1[1] = what_ref[1];
  eta_jp1[2] = what_ref[2];
  emlrtMEXProfilingStatement(15, isMexOutdated);
  what_j[0] = what_ref[0];
  what_j[1] = what_ref[1];
  what_j[2] = what_ref[2];
  /*  Calculate step-size (sigma from power iteration) */
  emlrtMEXProfilingStatement(16, isMexOutdated);
  st.site = &c_emlrtRSI;
  b_st.site = &q_emlrtRSI;
  c_st.site = &r_emlrtRSI;
  st.site = &c_emlrtRSI;
  t2 = lambda * lambda + 4.0 * omega * sigma;
  if (t2 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  t2 = muDoubleScalarSqrt(t2);
  alpha = 2.0 / (lambda + t2);
  emlrtMEXProfilingStatement(17, isMexOutdated);
  beta = omega * alpha;
  /*  Initialize */
  emlrtMEXProfilingStatement(18, isMexOutdated);
  memcpy(&zhat_j[0], &zeta_jp1[0], 1000U * sizeof(real_T));
  emlrtMEXProfilingStatement(19, isMexOutdated);
  memset(&zhat_jp1[0], 0, 1000U * sizeof(real_T));
  emlrtMEXProfilingStatement(20, isMexOutdated);
  what_star[0] = 0.0;
  what_star[1] = 0.0;
  what_star[2] = 0.0;
  emlrtMEXProfilingStatement(21, isMexOutdated);
  if ((int32_T)j_max != i) {
    emlrtIntegerCheckR2012b(j_max, &e_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = sol_info->zhis->size[0] * sol_info->zhis->size[1] *
      sol_info->zhis->size[2];
  sol_info->zhis->size[0] = 1000;
  sol_info->zhis->size[1] = 1;
  sol_info->zhis->size[2] = loop_ub_tmp;
  emxEnsureCapacity_real_T(sp, sol_info->zhis, i, &g_emlrtRTEI);
  loop_ub_tmp = 1000 * (int32_T)j_max;
  for (i = 0; i < loop_ub_tmp; i++) {
    sol_info->zhis->data[i] = 0.0;
  }
  emlrtMEXProfilingStatement(22, isMexOutdated);
  emlrtForLoopVectorCheckR2021a(1.0, 1.0, j_max, mxDOUBLE_CLASS, (int32_T)j_max,
                                &b_emlrtRTEI, (emlrtConstCTX)sp);
  loop_ub_tmp = 0;
  exitg1 = false;
  while ((!exitg1) && (loop_ub_tmp <= (int32_T)j_max - 1)) {
    __m128d r;
    __m128d r1;
    real_T dv[2];
    int32_T i2;
    int32_T i3;
    int32_T i4;
    int32_T i5;
    boolean_T guard1;
    /*     %% Projected gradient step */
    emlrtMEXProfilingStatement(23, isMexOutdated);
    st.site = &d_emlrtRSI;
    b_st.site = &n_emlrtRSI;
    b_mtimes(Hhat, eta_jp1, zhat_jp1);
    st.site = &d_emlrtRSI;
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zeta_jp1[i]);
      r1 = _mm_loadu_pd(&zhat_jp1[i]);
      _mm_storeu_pd(
          &b_zhat_jp1[i],
          _mm_sub_pd(
              r, _mm_mul_pd(
                     _mm_set1_pd(alpha),
                     _mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_set1_pd(lambda), r),
                                           _mm_loadu_pd(&qhat[i])),
                                r1))));
    }
    b_st.site = &n_emlrtRSI;
    mtimes(L_inv, b_zhat_jp1, zhat_jp1);
    st.site = &d_emlrtRSI;
    emlrtMEXProfilingFunctionEntry((char_T *)Ballz_project_complete_name,
                                   isMexOutdated);
    /* PROJECT Summary of this method goes here */
    /*    Detailed explanation goes here */
    emlrtMEXProfilingStatement(1, isMexOutdated);
    if (Dball->r <= 0.0) {
      b_st.site = &u_emlrtRSI;
      b_error(&b_st, e_emlrt_marshallOut(&b_st), &emlrtMCI);
    }
    emlrtMEXProfilingStatement(2, isMexOutdated);
    emlrtMEXProfilingStatement(5, isMexOutdated);
    emlrtMEXProfilingStatement(6, isMexOutdated);
    i = (int32_T)muDoubleScalarFloor(Dball->indices[0]);
    if (Dball->indices[0] != i) {
      emlrtIntegerCheckR2012b(Dball->indices[0], &l_emlrtDCI, &st);
    }
    if ((Dball->indices[0] < 1.0) || (Dball->indices[0] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dball->indices[0], 1, 1000,
                                    &i_emlrtBCI, &st);
    }
    i1 = (int32_T)muDoubleScalarFloor(Dball->indices[1]);
    if (Dball->indices[1] != i1) {
      emlrtIntegerCheckR2012b(Dball->indices[1], &l_emlrtDCI, &st);
    }
    if ((Dball->indices[1] < 1.0) || (Dball->indices[1] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dball->indices[1], 1, 1000,
                                    &i_emlrtBCI, &st);
    }
    emlrtMEXProfilingStatement(7, isMexOutdated);
    if ((int32_T)Dball->indices[0] != i) {
      emlrtIntegerCheckR2012b(Dball->indices[0], &k_emlrtDCI, &st);
    }
    if (((int16_T)Dball->indices[0] < 1) ||
        ((int16_T)Dball->indices[0] > 1000)) {
      emlrtDynamicBoundsCheckR2012b((int16_T)Dball->indices[0], 1, 1000,
                                    &h_emlrtBCI, &st);
    }
    if ((int32_T)Dball->indices[1] != i1) {
      emlrtIntegerCheckR2012b(Dball->indices[1], &k_emlrtDCI, &st);
    }
    if (((int16_T)Dball->indices[1] < 1) ||
        ((int16_T)Dball->indices[1] > 1000)) {
      emlrtDynamicBoundsCheckR2012b((int16_T)Dball->indices[1], 1, 1000,
                                    &h_emlrtBCI, &st);
    }
    dv[0] = zhat_jp1[(int16_T)Dball->indices[0] - 1];
    dv[1] = zhat_jp1[(int16_T)Dball->indices[1] - 1];
    emlrtMEXProfilingStatement(8, isMexOutdated);
    t2 = Dball->r / muDoubleScalarMax(vecnorm(dv), Dball->r);
    zhat_jp1[(int16_T)Dball->indices[0] - 1] = t2 * dv[0];
    zhat_jp1[(int16_T)Dball->indices[1] - 1] = t2 * dv[1];
    emlrtMEXProfilingStatement(9, isMexOutdated);
    emlrtMEXProfilingStatement(10, isMexOutdated);
    emlrtMEXProfilingFunctionExit(isMexOutdated);
    emlrtMEXProfilingStatement(24, isMexOutdated);
    st.site = &e_emlrtRSI;
    b_st.site = &e_emlrtRSI;
    emlrtMEXProfilingFunctionEntry((char_T *)c_Singletonz_project_complete_n,
                                   isMexOutdated);
    /* PROJECT Summary of this method goes here */
    /*    Detailed explanation goes here */
    emlrtMEXProfilingStatement(1, isMexOutdated);
    emlrtMEXProfilingStatement(4, isMexOutdated);
    emlrtMEXProfilingStatement(5, isMexOutdated);
    i = (int32_T)muDoubleScalarFloor(Dsingleton->indices[0]);
    if (Dsingleton->indices[0] != i) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[0], &j_emlrtDCI, &b_st);
    }
    if ((Dsingleton->indices[0] < 1.0) || (Dsingleton->indices[0] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[0], 1, 1000,
                                    &g_emlrtBCI, &b_st);
    }
    i1 = (int32_T)muDoubleScalarFloor(Dsingleton->indices[1]);
    if (Dsingleton->indices[1] != i1) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[1], &j_emlrtDCI, &b_st);
    }
    if ((Dsingleton->indices[1] < 1.0) || (Dsingleton->indices[1] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[1], 1, 1000,
                                    &g_emlrtBCI, &b_st);
    }
    i2 = (int32_T)muDoubleScalarFloor(Dsingleton->indices[2]);
    if (Dsingleton->indices[2] != i2) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[2], &j_emlrtDCI, &b_st);
    }
    if ((Dsingleton->indices[2] < 1.0) || (Dsingleton->indices[2] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[2], 1, 1000,
                                    &g_emlrtBCI, &b_st);
    }
    emlrtMEXProfilingStatement(6, isMexOutdated);
    i3 = (int32_T)Dsingleton->indices[0];
    if (i3 != i) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[0], &i_emlrtDCI, &b_st);
    }
    if ((Dsingleton->indices[0] < 1.0) || (Dsingleton->indices[0] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[0], 1, 1000,
                                    &f_emlrtBCI, &b_st);
    }
    i4 = (int32_T)Dsingleton->indices[1];
    if (i4 != i1) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[1], &i_emlrtDCI, &b_st);
    }
    if ((Dsingleton->indices[1] < 1.0) || (Dsingleton->indices[1] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[1], 1, 1000,
                                    &f_emlrtBCI, &b_st);
    }
    i5 = (int32_T)Dsingleton->indices[2];
    if (i5 != i2) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[2], &i_emlrtDCI, &b_st);
    }
    if ((Dsingleton->indices[2] < 1.0) || (Dsingleton->indices[2] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[2], 1, 1000,
                                    &f_emlrtBCI, &b_st);
    }
    emlrtMEXProfilingStatement(7, isMexOutdated);
    if (i3 != i) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[0], &h_emlrtDCI, &b_st);
    }
    if ((Dsingleton->indices[0] < 1.0) || (Dsingleton->indices[0] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[0], 1, 1000,
                                    &e_emlrtBCI, &b_st);
    }
    zhat_jp1[i3 - 1] = Dsingleton->y[0];
    if (i4 != i1) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[1], &h_emlrtDCI, &b_st);
    }
    if ((Dsingleton->indices[1] < 1.0) || (Dsingleton->indices[1] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[1], 1, 1000,
                                    &e_emlrtBCI, &b_st);
    }
    zhat_jp1[i4 - 1] = Dsingleton->y[1];
    if (i5 != i2) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[2], &h_emlrtDCI, &b_st);
    }
    if ((Dsingleton->indices[2] < 1.0) || (Dsingleton->indices[2] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[2], 1, 1000,
                                    &e_emlrtBCI, &b_st);
    }
    zhat_jp1[i5 - 1] = Dsingleton->y[2];
    emlrtMEXProfilingStatement(8, isMexOutdated);
    emlrtMEXProfilingStatement(9, isMexOutdated);
    emlrtMEXProfilingFunctionExit(isMexOutdated);
    memcpy(&b_zhat_jp1[0], &zhat_jp1[0], 1000U * sizeof(real_T));
    b_st.site = &n_emlrtRSI;
    mtimes(L, b_zhat_jp1, zhat_jp1);
    /* aff_violation(j) = norm(Hhat * zhat_jp1 - h); */
    emlrtMEXProfilingStatement(25, isMexOutdated);
    if (loop_ub_tmp + 1 > sol_info->zhis->size[2]) {
      emlrtDynamicBoundsCheckR2012b(loop_ub_tmp + 1, 1, sol_info->zhis->size[2],
                                    &b_emlrtBCI, (emlrtConstCTX)sp);
    }
    for (i = 0; i < 1000; i++) {
      sol_info->zhis->data[i + 1000 * loop_ub_tmp] = zhat_j[i];
    }
    /*     %% Proportional-integral feedback of affine equality constraint
     * violation */
    emlrtMEXProfilingStatement(26, isMexOutdated);
    st.site = &f_emlrtRSI;
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zhat_jp1[i]);
      r1 = _mm_loadu_pd(&zeta_jp1[i]);
      _mm_storeu_pd(&b_zhat_jp1[i],
                    _mm_sub_pd(_mm_mul_pd(_mm_set1_pd(2.0), r), r1));
    }
    b_st.site = &n_emlrtRSI;
    c_mtimes(Hhat, b_zhat_jp1, what_star);
    r = _mm_loadu_pd(&what_star[0]);
    r1 = _mm_loadu_pd(&eta_jp1[0]);
    _mm_storeu_pd(
        &what_star[0],
        _mm_add_pd(r1, _mm_mul_pd(_mm_set1_pd(beta),
                                  _mm_sub_pd(r, _mm_loadu_pd(&h[0])))));
    what_star[2] = eta_jp1[2] + beta * (what_star[2] - h[2]);
    /*     %% Extrapolate transformed primal variables */
    emlrtMEXProfilingStatement(27, isMexOutdated);
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zeta_jp1[i]);
      r1 = _mm_loadu_pd(&zhat_jp1[i]);
      _mm_storeu_pd(&zeta_jp1[i],
                    _mm_add_pd(_mm_mul_pd(_mm_set1_pd(1.0 - rho), r),
                               _mm_mul_pd(_mm_set1_pd(rho), r1)));
    }
    /*     %% Extrapolate transformed dual variable */
    emlrtMEXProfilingStatement(28, isMexOutdated);
    r = _mm_loadu_pd(&eta_jp1[0]);
    r1 = _mm_loadu_pd(&what_star[0]);
    _mm_storeu_pd(&eta_jp1[0], _mm_add_pd(_mm_mul_pd(_mm_set1_pd(1.0 - rho), r),
                                          _mm_mul_pd(_mm_set1_pd(rho), r1)));
    eta_jp1[2] = (1.0 - rho) * eta_jp1[2] + rho * what_star[2];
    /*     %% Check stopping criterion every j_check iterations */
    emlrtMEXProfilingStatement(29, isMexOutdated);
    guard1 = false;
    if (b_mod((real_T)loop_ub_tmp + 1.0, j_check) == 0.0) {
      real_T d;
      boolean_T terminate;
      emlrtMEXProfilingStatement(30, isMexOutdated);
      terminate = stopping(zhat_jp1, what_star, zhat_j, what_j, tol_abs,
                           tol_rel, &t2, &what_inf_dj);
      emlrtMEXProfilingStatement(31, isMexOutdated);
      d = ((real_T)loop_ub_tmp + 1.0) / j_check;
      i = (int32_T)muDoubleScalarFloor(d);
      if (d != i) {
        emlrtIntegerCheckR2012b(d, &f_emlrtDCI, (emlrtConstCTX)sp);
      }
      i1 = (int32_T)d;
      if ((d < 1.0) || (i1 > sol_info->zhat_inf_dj->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)d, 1,
                                      sol_info->zhat_inf_dj->size[0],
                                      &c_emlrtBCI, (emlrtConstCTX)sp);
      }
      sol_info->zhat_inf_dj->data[i1 - 1] = t2;
      emlrtMEXProfilingStatement(32, isMexOutdated);
      if (i1 != i) {
        emlrtIntegerCheckR2012b(d, &g_emlrtDCI, (emlrtConstCTX)sp);
      }
      if ((d < 1.0) || (i1 > sol_info->what_inf_dj->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)d, 1,
                                      sol_info->what_inf_dj->size[0],
                                      &d_emlrtBCI, (emlrtConstCTX)sp);
      }
      sol_info->what_inf_dj->data[i1 - 1] = what_inf_dj;
      emlrtMEXProfilingStatement(33, isMexOutdated);
      if (terminate) {
        /*  stopping(zhat_jp1, what_jp1, zhat_j, what_j, tol_abs, tol_rel) */
        emlrtMEXProfilingStatement(34, isMexOutdated);
        sol_info->iterations = (real_T)loop_ub_tmp + 1.0;
        emlrtMEXProfilingStatement(35, isMexOutdated);
        sol_info->solution_status.Value.size[0] = 1;
        sol_info->solution_status.Value.size[1] = 7;
        for (i = 0; i < 7; i++) {
          sol_info->solution_status.Value.data[i] = cv1[i];
        }
        /* fprintf("Opt at %d\n", j) */
        emlrtMEXProfilingStatement(36, isMexOutdated);
        exitg1 = true;
      } else {
        emlrtMEXProfilingStatement(37, isMexOutdated);
        /*  Try to check infeasibility - should be test which can be applied
         * before end... */
        emlrtMEXProfilingStatement(38, isMexOutdated);
        terminate = false;
        emlrtMEXProfilingStatement(39, isMexOutdated);
        if ((what_inf_dj < tol_abs) && (t2 > tol_infeas)) {
          emlrtMEXProfilingStatement(40, isMexOutdated);
          sol_info->primal_infeasible = true;
          emlrtMEXProfilingStatement(41, isMexOutdated);
          sol_info->iterations = (real_T)loop_ub_tmp + 1.0;
          emlrtMEXProfilingStatement(42, isMexOutdated);
          sol_info->solution_status.Value.size[0] = 1;
          sol_info->solution_status.Value.size[1] = 10;
          for (i = 0; i < 10; i++) {
            sol_info->solution_status.Value.data[i] = cv2[i];
          }
          emlrtMEXProfilingStatement(43, isMexOutdated);
          terminate = true;
          /* fprintf("Prim infeas at %d\n", j) */
          emlrtMEXProfilingStatement(44, isMexOutdated);
        }
        emlrtMEXProfilingStatement(45, isMexOutdated);
        if ((t2 < tol_abs) && (what_inf_dj > tol_infeas)) {
          emlrtMEXProfilingStatement(46, isMexOutdated);
          sol_info->dual_infeasible = true;
          emlrtMEXProfilingStatement(47, isMexOutdated);
          sol_info->iterations = (real_T)loop_ub_tmp + 1.0;
          emlrtMEXProfilingStatement(48, isMexOutdated);
          sol_info->solution_status.Value.size[0] = 1;
          sol_info->solution_status.Value.size[1] = 10;
          for (i = 0; i < 10; i++) {
            sol_info->solution_status.Value.data[i] = cv2[i];
          }
          emlrtMEXProfilingStatement(49, isMexOutdated);
          terminate = true;
          /* fprintf("Dual infeas at %d\n", j) */
          emlrtMEXProfilingStatement(50, isMexOutdated);
        }
        emlrtMEXProfilingStatement(51, isMexOutdated);
        if (terminate) {
          emlrtMEXProfilingStatement(52, isMexOutdated);
          exitg1 = true;
        } else {
          emlrtMEXProfilingStatement(54, isMexOutdated);
          emlrtMEXProfilingStatement(55, isMexOutdated);
          guard1 = true;
        }
      }
    } else {
      guard1 = true;
    }
    if (guard1) {
      /*     %% Update values for next iteration */
      emlrtMEXProfilingStatement(56, isMexOutdated);
      memcpy(&zhat_j[0], &zhat_jp1[0], 1000U * sizeof(real_T));
      emlrtMEXProfilingStatement(57, isMexOutdated);
      what_j[0] = what_star[0];
      what_j[1] = what_star[1];
      what_j[2] = what_star[2];
      emlrtMEXProfilingStatement(58, isMexOutdated);
      emlrtMEXProfilingStatement(59, isMexOutdated);
      emlrtMEXProfilingStatement(60, isMexOutdated);
      loop_ub_tmp++;
    }
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  /*  Recover original primal variables */
  emlrtMEXProfilingStatement(61, isMexOutdated);
  st.site = &g_emlrtRSI;
  b_st.site = &n_emlrtRSI;
  mtimes(L_inv, zhat_jp1, z_star);
  /*  Retain transformed dual variable */
  emlrtMEXProfilingStatement(62, isMexOutdated);
  emlrtMEXProfilingStatement(63, isMexOutdated);
  st.site = &h_emlrtRSI;
  t2 = toc(&st, expl_temp.tv_sec, expl_temp.tv_nsec);
  emlrtMEXProfilingStatement(64, isMexOutdated);
  sol_info->time = t2;
  /*  %  */
  emlrtMEXProfilingStatement(65, isMexOutdated);
  st.site = &i_emlrtRSI;
  b_st.site = &t_emlrtRSI;
  y = NULL;
  m = emlrtCreateCharArray(2, &iv[0]);
  emlrtInitCharArrayR2013a(&b_st, 7, m, &u[0]);
  emlrtAssign(&y, m);
  b_y = NULL;
  m = emlrtCreateDoubleScalar(1.0);
  emlrtAssign(&b_y, m);
  c_y = NULL;
  m = emlrtCreateCharArray(2, &iv1[0]);
  emlrtInitCharArrayR2013a(&b_st, 20, m, &b_u[0]);
  emlrtAssign(&c_y, m);
  d_y = NULL;
  m = emlrtCreateDoubleScalar(t2 * 1000.0);
  emlrtAssign(&d_y, m);
  c_st.site = &v_emlrtRSI;
  t2 = emlrt_marshallIn(&c_st, feval(&c_st, y, b_y, c_y, d_y, &b_emlrtMCI),
                        "<output of feval>");
  e_y = NULL;
  m = emlrtCreateDoubleScalar(t2);
  emlrtAssign(&e_y, m);
  emlrtDisplayR2012b(e_y, "ans", &emlrtRTEI, (emlrtCTX)sp);
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
  emlrtMEXProfilingStatement(66, isMexOutdated);
  emlrtMEXProfilingStatement(67, isMexOutdated);
  emlrtMEXProfilingFunctionExit(isMexOutdated);
}

/* End of code generation (PIPG_new_explicit.c) */
