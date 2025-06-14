/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_orig_precond_Dexplicit.c
 *
 * Code generation for function 'PIPG_orig_precond_Dexplicit'
 *
 */

/* Include files */
#include "PIPG_orig_precond_Dexplicit.h"
#include "PIPG_orig_precond_Dexplicit_data.h"
#include "PIPG_orig_precond_Dexplicit_emxutil.h"
#include "PIPG_orig_precond_Dexplicit_mexutil.h"
#include "PIPG_orig_precond_Dexplicit_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "mod.h"
#include "mtimes.h"
#include "mtimes1.h"
#include "rt_nonfinite.h"
#include "sparse.h"
#include "tic.h"
#include "toc.h"
#include "blas.h"
#include "emlrt.h"
#include "mwmathutil.h"
#include <emmintrin.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Variable Definitions */
static emlrtRTEInfo emlrtRTEI = {
    132,                           /* lineNo */
    1,                             /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pName */
};

static emlrtRSInfo emlrtRSI = {
    21,                            /* lineNo */
    "PIPG_orig_precond_Dexplicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    35,                            /* lineNo */
    "PIPG_orig_precond_Dexplicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    52,                            /* lineNo */
    "PIPG_orig_precond_Dexplicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    55,                            /* lineNo */
    "PIPG_orig_precond_Dexplicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    56,                            /* lineNo */
    "PIPG_orig_precond_Dexplicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    62,                            /* lineNo */
    "PIPG_orig_precond_Dexplicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    125,                           /* lineNo */
    "PIPG_orig_precond_Dexplicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    130,                           /* lineNo */
    "PIPG_orig_precond_Dexplicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI = {
    132,                           /* lineNo */
    "PIPG_orig_precond_Dexplicit", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pathName */
};

static emlrtRSInfo n_emlrtRSI =
    {
        94,                  /* lineNo */
        "eml_mtimes_helper", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_"
        "helper.m" /* pathName */
};

static emlrtRSInfo y_emlrtRSI = {
    18,              /* lineNo */
    "sparse/mtimes", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\mtimes.m" /* pathName */
};

static emlrtRSInfo ob_emlrtRSI = {
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
    17,                            /* lineNo */
    36,                            /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    1                                     /* checkKind */
};

static emlrtBCInfo emlrtBCI = {
    1,                             /* iFirst */
    1,                             /* iLast */
    17,                            /* lineNo */
    36,                            /* colNo */
    "[size(z_ref, j_max)]",        /* aName */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    0                                     /* checkKind */
};

static emlrtRTEInfo b_emlrtRTEI = {
    47,                            /* lineNo */
    9,                             /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pName */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,                            /* iFirst */
    -1,                            /* iLast */
    59,                            /* lineNo */
    13,                            /* colNo */
    "zhis",                        /* aName */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    0                                     /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = {
    15,                            /* lineNo */
    30,                            /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    4                                     /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = {
    15,                            /* lineNo */
    30,                            /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    1                                     /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = {
    16,                            /* lineNo */
    30,                            /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    1                                     /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = {
    45,                            /* lineNo */
    14,                            /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    1                                     /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,                            /* iFirst */
    -1,                            /* iLast */
    73,                            /* lineNo */
    30,                            /* colNo */
    "sol_info.zhat_inf_dj",        /* aName */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    0                                     /* checkKind */
};

static emlrtDCInfo f_emlrtDCI = {
    73,                            /* lineNo */
    30,                            /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    1                                     /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,                            /* iFirst */
    -1,                            /* iLast */
    74,                            /* lineNo */
    30,                            /* colNo */
    "sol_info.what_inf_dj",        /* aName */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    0                                     /* checkKind */
};

static emlrtDCInfo g_emlrtDCI = {
    74,                            /* lineNo */
    30,                            /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m", /* pName */
    1                                     /* checkKind */
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

static emlrtRTEInfo i_emlrtRTEI = {
    15,                            /* lineNo */
    1,                             /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pName */
};

static emlrtRTEInfo j_emlrtRTEI = {
    16,                            /* lineNo */
    1,                             /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pName */
};

static emlrtRTEInfo k_emlrtRTEI = {
    45,                            /* lineNo */
    1,                             /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pName */
};

static emlrtRTEInfo l_emlrtRTEI = {
    55,                            /* lineNo */
    98,                            /* colNo */
    "PIPG_orig_precond_Dexplicit", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\PIPG_orig_precond_Dexplicit.m" /* pName */
};

static emlrtRSInfo pb_emlrtRSI = {
    23,              /* lineNo */
    "Ballz/project", /* fcnName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
    "IPG\\Convex Cones\\Ballz.m" /* pathName */
};

static emlrtRSInfo qb_emlrtRSI = {
    66,        /* lineNo */
    "fprintf", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\iofun\\fprintf.m" /* pathName
                                                                          */
};

/* Function Declarations */
static void b_error(const emlrtStack *sp, const mxArray *m,
                    emlrtMCInfo *location);

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

void PIPG_orig_precond_Dexplicit(
    c_PIPG_orig_precond_DexplicitSt *SD, const emlrtStack *sp,
    const real_T qhat[1000], const sparse Hhat, const real_T h[3],
    const Ballz *Dball, const Singletonz *Dsingleton, const real_T L[1000000],
    const real_T L_inv[1000000], real_T lambda, real_T sigma, real_T rho,
    real_T tol_abs, real_T tol_rel, real_T tol_infeas, real_T j_check,
    real_T j_restart, real_T j_max, const real_T z_ref[1000],
    const real_T what_ref[3], real_T S_z, real_T c_z, real_T z_star[1000],
    real_T what_star[3], struct0_T *sol_info)
{
  static const int32_T iv[2] = {1, 7};
  static const int32_T iv1[2] = {1, 20};
  static const int32_T iv2[2] = {1, 25};
  static const char_T c_u[25] = {'R', 'a', 'd', 'i', 'u', 's', ' ', 's', 'h',
                                 'o', 'u', 'l', 'd', ' ', 'b', 'e', ' ', 'p',
                                 'o', 's', 'i', 't', 'i', 'v', 'e'};
  static const char_T b_u[20] = {'P', 'I', 'P', 'G', ' ',  'T', 'i',
                                 'm', 'e', ':', ' ', '%',  '.', '3',
                                 'f', ' ', 'm', 's', '\\', 'n'};
  static const char_T cv2[10] = {'I', 'n', 'f', 'e', 'a',
                                 's', 'i', 'b', 'l', 'e'};
  static const char_T b_cv[8] = {'U', 'n', 's', 'o', 'l', 'v', 'e', 'd'};
  static const char_T cv1[7] = {'O', 'p', 't', 'i', 'm', 'a', 'l'};
  static const char_T u[7] = {'f', 'p', 'r', 'i', 'n', 't', 'f'};
  __m128d r;
  ptrdiff_t incx_t;
  ptrdiff_t incy_t;
  ptrdiff_t n_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  emlrtTimespec expl_temp;
  emxArray_int32_T *a_colidx;
  emxArray_int32_T *a_rowidx;
  emxArray_real_T *a_d;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *m;
  const mxArray *propValues;
  const mxArray *y;
  real_T b_zhat_j[1000];
  real_T c[1000];
  real_T zhat_j[1000];
  real_T zhat_jp1[1000];
  real_T vhat_j[3];
  real_T vhat_jp1[3];
  real_T what_j[3];
  real_T bc;
  real_T j_r;
  real_T scale;
  real_T *a_d_data;
  int32_T acol;
  int32_T ap;
  int32_T i;
  int32_T i1;
  int32_T j;
  int32_T nap;
  int32_T *a_colidx_data;
  int32_T *a_rowidx_data;
  const char_T *propClasses = "coder.internal.string";
  const char_T *propNames = "Value";
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /* PIPG Proportional Integral Projected Gradient convex conic solver */
  /*    Extrapolated Proportional Integral Projected Gradient (xPIPG) */
  /*  rho in [1.5, 1.9] is usually a good choice - should converge about twice
   */
  /*  as fast as not doing extrapolation (rho = 1) */
  /*  sigma should be >= norm(H) */
  /*  lambda should be >= norm(Q) */
  sol_info->primal_infeasible = false;
  sol_info->dual_infeasible = false;
  sol_info->solution_status.Value.size[0] = 1;
  sol_info->solution_status.Value.size[1] = 8;
  for (i = 0; i < 8; i++) {
    sol_info->solution_status.Value.data[i] = b_cv[i];
  }
  sol_info->iterations = -1.0;
  if (!(j_max >= 0.0)) {
    emlrtNonNegativeCheckR2012b(j_max, &b_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = (int32_T)muDoubleScalarFloor(j_max);
  if (j_max != i) {
    emlrtIntegerCheckR2012b(j_max, &c_emlrtDCI, (emlrtConstCTX)sp);
  }
  nap = (int32_T)j_max;
  i1 = sol_info->zhat_inf_dj->size[0];
  sol_info->zhat_inf_dj->size[0] = nap;
  emxEnsureCapacity_real_T(sp, sol_info->zhat_inf_dj, i1, &i_emlrtRTEI);
  for (i1 = 0; i1 < nap; i1++) {
    sol_info->zhat_inf_dj->data[i1] = 0.0;
  }
  if (nap != i) {
    emlrtIntegerCheckR2012b(j_max, &d_emlrtDCI, (emlrtConstCTX)sp);
  }
  i1 = sol_info->what_inf_dj->size[0];
  sol_info->what_inf_dj->size[0] = nap;
  emxEnsureCapacity_real_T(sp, sol_info->what_inf_dj, i1, &j_emlrtRTEI);
  for (i1 = 0; i1 < nap; i1++) {
    sol_info->what_inf_dj->data[i1] = 0.0;
  }
  if (nap != i) {
    emlrtIntegerCheckR2012b(j_max, &emlrtDCI, (emlrtConstCTX)sp);
  }
  if (j_max < 1.0) {
    emlrtDynamicBoundsCheckR2012b(0, 1, 1, &emlrtBCI, (emlrtConstCTX)sp);
  }
  st.site = &emlrtRSI;
  expl_temp = tic(&st);
  /*  Intialize variables */
  r = _mm_mul_pd(_mm_loadu_pd(&what_ref[0]), _mm_set1_pd(0.0));
  _mm_storeu_pd(&vhat_j[0], r);
  _mm_storeu_pd(&what_j[0], r);
  bc = what_ref[2] * 0.0;
  vhat_j[2] = bc;
  what_j[2] = bc;
  /*  Calculate step-size (sigma from power iteration) */
  /* alpha = 2 / (lambda + sqrt(lambda^2 + 4 * sigma)); */
  /*  Initialize */
  st.site = &b_emlrtRSI;
  b_st.site = &n_emlrtRSI;
  mtimes(L, z_ref, zhat_j);
  memset(&zhat_jp1[0], 0, 1000U * sizeof(real_T));
  what_star[0] = 0.0;
  what_star[1] = 0.0;
  what_star[2] = 0.0;
  j_r = 1.0;
  /*   */
  /*  zhat_inf_dj2 = 1e5; */
  /*  what_inf_dj2 = 1e5; */
  if ((int32_T)j_max != i) {
    emlrtIntegerCheckR2012b(j_max, &e_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = sol_info->zhis->size[0] * sol_info->zhis->size[1] *
      sol_info->zhis->size[2];
  sol_info->zhis->size[0] = 1000;
  sol_info->zhis->size[1] = 1;
  sol_info->zhis->size[2] = (int32_T)j_max;
  emxEnsureCapacity_real_T(sp, sol_info->zhis, i, &k_emlrtRTEI);
  nap = 1000 * (int32_T)j_max;
  for (i = 0; i < nap; i++) {
    sol_info->zhis->data[i] = 0.0;
  }
  emlrtForLoopVectorCheckR2021a(1.0, 1.0, j_max, mxDOUBLE_CLASS, (int32_T)j_max,
                                &b_emlrtRTEI, (emlrtConstCTX)sp);
  j = 0;
  emxInit_real_T(sp, &a_d, 1, &l_emlrtRTEI);
  emxInit_int32_T(sp, &a_colidx, &l_emlrtRTEI);
  emxInit_int32_T(sp, &a_rowidx, &l_emlrtRTEI);
  exitg1 = false;
  while ((!exitg1) && (j <= (int32_T)j_max - 1)) {
    __m128d r1;
    __m128d r2;
    __m128d r3;
    real_T absx;
    real_T absxk_tmp;
    real_T beta_j;
    real_T t;
    real_T varargin_1;
    int32_T apend1;
    boolean_T guard1;
    scale = 2.0 / (2.0 * lambda + (j_r + 1.0));
    beta_j = (j_r + 1.0) / (2.0 * sigma);
    /*     %% Proportional-integral feedback of affine equality constraint
     * violation */
    st.site = &c_emlrtRSI;
    sparse_mtimes(&st, Hhat.d, Hhat.colidx, Hhat.rowidx, zhat_j, what_star);
    r = _mm_loadu_pd(&what_star[0]);
    r1 = _mm_loadu_pd(&vhat_j[0]);
    _mm_storeu_pd(
        &what_star[0],
        _mm_add_pd(r1, _mm_mul_pd(_mm_set1_pd(beta_j),
                                  _mm_sub_pd(r, _mm_loadu_pd(&h[0])))));
    what_star[2] = vhat_j[2] + beta_j * (what_star[2] - h[2]);
    /*     %% Projected gradient step */
    st.site = &d_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    sparse_ctranspose(&b_st, Hhat.d, Hhat.colidx, Hhat.rowidx, a_d, a_colidx,
                      a_rowidx);
    a_rowidx_data = a_rowidx->data;
    a_colidx_data = a_colidx->data;
    a_d_data = a_d->data;
    b_st.site = &q_emlrtRSI;
    memset(&c[0], 0, 1000U * sizeof(real_T));
    if (a_colidx_data[a_colidx->size[0] - 1] - 1 != 0) {
      c_st.site = &r_emlrtRSI;
      memset(&c[0], 0, 1000U * sizeof(real_T));
      for (acol = 0; acol < 3; acol++) {
        bc = what_star[acol];
        i = a_colidx_data[acol];
        i1 = a_colidx_data[acol + 1];
        nap = i1 - a_colidx_data[acol];
        if (nap >= 4) {
          apend1 = ((i1 - nap) + ((nap / 4) << 2)) - 1;
          d_st.site = &s_emlrtRSI;
          if ((a_colidx_data[acol] <= apend1) && (apend1 > 2147483643)) {
            e_st.site = &t_emlrtRSI;
            check_forloop_overflow_error(&e_st);
          }
          for (ap = i; ap <= apend1; ap += 4) {
            nap = a_rowidx_data[ap - 1] - 1;
            c[nap] += a_d_data[ap - 1] * bc;
            c[a_rowidx_data[ap] - 1] += a_d_data[ap] * bc;
            nap = a_rowidx_data[ap + 1] - 1;
            c[nap] += a_d_data[ap + 1] * bc;
            nap = a_rowidx_data[ap + 2] - 1;
            c[nap] += a_d_data[ap + 2] * bc;
          }
          i = apend1 + 1;
          for (ap = i; ap < i1; ap++) {
            nap = a_rowidx_data[ap - 1] - 1;
            c[nap] += a_d_data[ap - 1] * bc;
          }
        } else {
          for (ap = i; ap < i1; ap++) {
            nap = a_rowidx_data[ap - 1] - 1;
            c[nap] += a_d_data[ap - 1] * bc;
          }
        }
      }
    }
    st.site = &d_emlrtRSI;
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zhat_j[i]);
      r1 = _mm_loadu_pd(&c[i]);
      _mm_storeu_pd(
          &b_zhat_j[i],
          _mm_sub_pd(
              r, _mm_mul_pd(
                     _mm_set1_pd(scale),
                     _mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_set1_pd(lambda), r),
                                           _mm_loadu_pd(&qhat[i])),
                                r1))));
    }
    b_st.site = &n_emlrtRSI;
    mtimes(L_inv, b_zhat_j, zhat_jp1);
    st.site = &d_emlrtRSI;
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zhat_jp1[i]);
      _mm_storeu_pd(&zhat_jp1[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd(S_z), r),
                                             _mm_set1_pd(c_z)));
    }
    /* PROJECT Summary of this method goes here */
    /*    Detailed explanation goes here */
    if (Dball->r <= 0.0) {
      f_y = NULL;
      m = NULL;
      emlrtAssign(&f_y,
                  emlrtCreateClassInstance2022a(&st, "coder.internal.string"));
      m = NULL;
      g_y = NULL;
      propValues = emlrtCreateCharArray(2, &iv2[0]);
      emlrtInitCharArrayR2013a(&st, 25, propValues, &c_u[0]);
      emlrtAssign(&g_y, propValues);
      emlrtAssign(&m, g_y);
      propValues = m;
      emlrtSetAllProperties(&st, &f_y, 0, 1, (const char_T **)&propNames,
                            (const char_T **)&propClasses, &propValues);
      emlrtAssign(&f_y, emlrtConvertInstanceToRedirectSource(
                            &st, f_y, 0, "coder.internal.string"));
      b_st.site = &pb_emlrtRSI;
      b_error(&b_st, f_y, &emlrtMCI);
    }
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
    scale = 3.3121686421112381E-170;
    absx = zhat_jp1[(int16_T)Dball->indices[0] - 1];
    bc = muDoubleScalarAbs(absx);
    if (bc > 3.3121686421112381E-170) {
      varargin_1 = 1.0;
      scale = bc;
    } else {
      t = bc / 3.3121686421112381E-170;
      varargin_1 = t * t;
    }
    absxk_tmp = zhat_jp1[(int16_T)Dball->indices[1] - 1];
    bc = muDoubleScalarAbs(absxk_tmp);
    if (bc > scale) {
      t = scale / bc;
      varargin_1 = varargin_1 * t * t + 1.0;
      scale = bc;
    } else {
      t = bc / scale;
      varargin_1 += t * t;
    }
    scale = Dball->r /
            muDoubleScalarMax(scale * muDoubleScalarSqrt(varargin_1), Dball->r);
    zhat_jp1[(int16_T)Dball->indices[0] - 1] = scale * absx;
    zhat_jp1[(int16_T)Dball->indices[1] - 1] = scale * absxk_tmp;
    st.site = &e_emlrtRSI;
    /* PROJECT Summary of this method goes here */
    /*    Detailed explanation goes here */
    i = (int32_T)muDoubleScalarFloor(Dsingleton->indices[0]);
    if (Dsingleton->indices[0] != i) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[0], &j_emlrtDCI, &st);
    }
    if ((Dsingleton->indices[0] < 1.0) || (Dsingleton->indices[0] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[0], 1, 1000,
                                    &g_emlrtBCI, &st);
    }
    i1 = (int32_T)muDoubleScalarFloor(Dsingleton->indices[1]);
    if (Dsingleton->indices[1] != i1) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[1], &j_emlrtDCI, &st);
    }
    if ((Dsingleton->indices[1] < 1.0) || (Dsingleton->indices[1] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[1], 1, 1000,
                                    &g_emlrtBCI, &st);
    }
    nap = (int32_T)muDoubleScalarFloor(Dsingleton->indices[2]);
    if (Dsingleton->indices[2] != nap) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[2], &j_emlrtDCI, &st);
    }
    if ((Dsingleton->indices[2] < 1.0) || (Dsingleton->indices[2] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[2], 1, 1000,
                                    &g_emlrtBCI, &st);
    }
    apend1 = (int32_T)Dsingleton->indices[0];
    if (apend1 != i) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[0], &i_emlrtDCI, &st);
    }
    if ((Dsingleton->indices[0] < 1.0) || (Dsingleton->indices[0] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[0], 1, 1000,
                                    &f_emlrtBCI, &st);
    }
    ap = (int32_T)Dsingleton->indices[1];
    if (ap != i1) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[1], &i_emlrtDCI, &st);
    }
    if ((Dsingleton->indices[1] < 1.0) || (Dsingleton->indices[1] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[1], 1, 1000,
                                    &f_emlrtBCI, &st);
    }
    acol = (int32_T)Dsingleton->indices[2];
    if (acol != nap) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[2], &i_emlrtDCI, &st);
    }
    if ((Dsingleton->indices[2] < 1.0) || (Dsingleton->indices[2] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[2], 1, 1000,
                                    &f_emlrtBCI, &st);
    }
    if (apend1 != i) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[0], &h_emlrtDCI, &st);
    }
    if ((Dsingleton->indices[0] < 1.0) || (Dsingleton->indices[0] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[0], 1, 1000,
                                    &e_emlrtBCI, &st);
    }
    zhat_jp1[apend1 - 1] = Dsingleton->y[0];
    if (ap != i1) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[1], &h_emlrtDCI, &st);
    }
    if ((Dsingleton->indices[1] < 1.0) || (Dsingleton->indices[1] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[1], 1, 1000,
                                    &e_emlrtBCI, &st);
    }
    zhat_jp1[ap - 1] = Dsingleton->y[1];
    if (acol != nap) {
      emlrtIntegerCheckR2012b(Dsingleton->indices[2], &h_emlrtDCI, &st);
    }
    if ((Dsingleton->indices[2] < 1.0) || (Dsingleton->indices[2] > 1000.0)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)Dsingleton->indices[2], 1, 1000,
                                    &e_emlrtBCI, &st);
    }
    zhat_jp1[acol - 1] = Dsingleton->y[2];
    st.site = &e_emlrtRSI;
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zhat_jp1[i]);
      _mm_storeu_pd(&b_zhat_j[i], _mm_div_pd(_mm_sub_pd(r, _mm_set1_pd(c_z)),
                                             _mm_set1_pd(S_z)));
    }
    b_st.site = &n_emlrtRSI;
    mtimes(L, b_zhat_j, zhat_jp1);
    /*   */
    /*  aff_violation(j) = norm(Hhat * zhat_jp1 - h); */
    if (j + 1 > sol_info->zhis->size[2]) {
      emlrtDynamicBoundsCheckR2012b(j + 1, 1, sol_info->zhis->size[2],
                                    &b_emlrtBCI, (emlrtConstCTX)sp);
    }
    for (i = 0; i < 1000; i++) {
      sol_info->zhis->data[i + 1000 * j] = zhat_j[i];
    }
    /*     %% */
    st.site = &f_emlrtRSI;
    b_st.site = &y_emlrtRSI;
    sparse_times(&b_st, beta_j, Hhat.d, Hhat.colidx, Hhat.rowidx, a_d, a_colidx,
                 a_rowidx);
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zhat_jp1[i]);
      r1 = _mm_loadu_pd(&zhat_j[i]);
      _mm_storeu_pd(&b_zhat_j[i], _mm_sub_pd(r, r1));
    }
    st.site = &f_emlrtRSI;
    sparse_mtimes(&st, a_d, a_colidx, a_rowidx, b_zhat_j, vhat_jp1);
    /*     %% Extrapolate variables */
    r = _mm_loadu_pd(&what_star[0]);
    r1 = _mm_loadu_pd(&vhat_jp1[0]);
    _mm_storeu_pd(&vhat_jp1[0], _mm_add_pd(r, r1));
    r1 = _mm_loadu_pd(&what_j[0]);
    r2 = _mm_set1_pd(1.0 - rho);
    r3 = _mm_set1_pd(rho);
    _mm_storeu_pd(&what_star[0],
                  _mm_add_pd(_mm_mul_pd(r2, r1), _mm_mul_pd(r3, r)));
    vhat_jp1[2] += what_star[2];
    what_star[2] = (1.0 - rho) * what_j[2] + rho * what_star[2];
    for (i = 0; i <= 998; i += 2) {
      r = _mm_loadu_pd(&zhat_j[i]);
      r1 = _mm_loadu_pd(&zhat_jp1[i]);
      _mm_storeu_pd(&zhat_jp1[i],
                    _mm_add_pd(_mm_mul_pd(r2, r), _mm_mul_pd(r3, r1)));
    }
    r = _mm_loadu_pd(&vhat_j[0]);
    r1 = _mm_loadu_pd(&vhat_jp1[0]);
    _mm_storeu_pd(&vhat_jp1[0],
                  _mm_add_pd(_mm_mul_pd(r2, r), _mm_mul_pd(r3, r1)));
    vhat_jp1[2] = (1.0 - rho) * vhat_j[2] + rho * vhat_jp1[2];
    /*     %% Check stopping criterion every j_check iterations */
    guard1 = false;
    if (b_mod((real_T)j + 1.0, j_check) == 0.0) {
      real_T what_inf_dj;
      real_T what_inf_j;
      boolean_T terminate;
      /* STOPPING Summary of this function goes here */
      /*    Detailed explanation goes here */
      scale = 0.0;
      t = 0.0;
      varargin_1 = 0.0;
      for (nap = 0; nap < 1000; nap++) {
        bc = zhat_jp1[nap];
        absx = muDoubleScalarAbs(bc);
        if (muDoubleScalarIsNaN(absx) || (absx > scale)) {
          scale = absx;
        }
        absxk_tmp = zhat_j[nap];
        absx = muDoubleScalarAbs(absxk_tmp);
        if (muDoubleScalarIsNaN(absx) || (absx > t)) {
          t = absx;
        }
        bc -= absxk_tmp;
        c[nap] = bc;
        absx = muDoubleScalarAbs(bc);
        if (muDoubleScalarIsNaN(absx) || (absx > varargin_1)) {
          varargin_1 = absx;
        }
      }
      beta_j = 0.0;
      what_inf_j = 0.0;
      what_inf_dj = 0.0;
      for (nap = 0; nap < 3; nap++) {
        bc = what_star[nap];
        absx = muDoubleScalarAbs(bc);
        if (muDoubleScalarIsNaN(absx) || (absx > beta_j)) {
          beta_j = absx;
        }
        absxk_tmp = what_j[nap];
        absx = muDoubleScalarAbs(absxk_tmp);
        if (muDoubleScalarIsNaN(absx) || (absx > what_inf_j)) {
          what_inf_j = absx;
        }
        bc -= absxk_tmp;
        what_j[nap] = bc;
        absx = muDoubleScalarAbs(bc);
        if (muDoubleScalarIsNaN(absx) || (absx > what_inf_dj)) {
          what_inf_dj = absx;
        }
      }
      if ((varargin_1 <= tol_abs + tol_rel * muDoubleScalarMax(scale, t)) &&
          (what_inf_dj <=
           tol_abs + tol_rel * muDoubleScalarMax(beta_j, what_inf_j))) {
        terminate = true;
      } else {
        terminate = false;
      }
      bc = ((real_T)j + 1.0) / j_check;
      i = (int32_T)muDoubleScalarFloor(bc);
      if (bc != i) {
        emlrtIntegerCheckR2012b(bc, &f_emlrtDCI, (emlrtConstCTX)sp);
      }
      i1 = (int32_T)bc;
      if ((bc < 1.0) || (i1 > sol_info->zhat_inf_dj->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)bc, 1,
                                      sol_info->zhat_inf_dj->size[0],
                                      &c_emlrtBCI, (emlrtConstCTX)sp);
      }
      sol_info->zhat_inf_dj->data[i1 - 1] = varargin_1;
      if (i1 != i) {
        emlrtIntegerCheckR2012b(bc, &g_emlrtDCI, (emlrtConstCTX)sp);
      }
      if ((bc < 1.0) || (i1 > sol_info->what_inf_dj->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)bc, 1,
                                      sol_info->what_inf_dj->size[0],
                                      &d_emlrtBCI, (emlrtConstCTX)sp);
      }
      sol_info->what_inf_dj->data[i1 - 1] = what_inf_dj;
      if (terminate) {
        /*  stopping(zhat_jp1, what_jp1, zhat_j, what_j, tol_abs, tol_rel) */
        sol_info->iterations = (real_T)j + 1.0;
        sol_info->solution_status.Value.size[0] = 1;
        sol_info->solution_status.Value.size[1] = 7;
        for (i = 0; i < 7; i++) {
          sol_info->solution_status.Value.data[i] = cv1[i];
        }
        /* fprintf("Opt at %d\n", j) */
        exitg1 = true;
      } else {
        /*  Try to check infeasibility - should be test which can be applied
         * before the end... */
        if ((what_inf_dj < tol_abs) && (varargin_1 > tol_infeas)) {
          sol_info->primal_infeasible = true;
          sol_info->iterations = (real_T)j + 1.0;
          sol_info->solution_status.Value.size[0] = 1;
          sol_info->solution_status.Value.size[1] = 10;
          for (i = 0; i < 10; i++) {
            sol_info->solution_status.Value.data[i] = cv2[i];
          }
          /* terminate = true; */
          /* fprintf("Prim infeas at %d\n", j) */
        }
        if ((varargin_1 < tol_abs) && (what_inf_dj > tol_infeas)) {
          sol_info->dual_infeasible = true;
          sol_info->iterations = (real_T)j + 1.0;
          sol_info->solution_status.Value.size[0] = 1;
          sol_info->solution_status.Value.size[1] = 10;
          for (i = 0; i < 10; i++) {
            sol_info->solution_status.Value.data[i] = cv2[i];
          }
          /* terminate = true; */
          /* fprintf("Dual infeas at %d\n", j) */
        }
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
    if (guard1) {
      /*     %% Restart */
      /*  [~, zhat_inf_djp1, what_inf_djp1] = stopping(zhat_jp1, what_jp1,
       * zhat_j, what_j, tol_abs, tol_rel); */
      /*     if (zhat_inf_djp1 - zhat_inf_dj2  <= 0) %(vhat_jp1 - zhat_jp1)' *
       * (zhat_jp1 - zhat_j) <= 0 && j_restart == -1 || mod(j, j_restart) ~= 0
       */
      r = _mm_loadu_pd(&vhat_jp1[0]);
      r1 = _mm_loadu_pd(&zhat_jp1[0]);
      _mm_storeu_pd(&c[0], _mm_sub_pd(r, r1));
      c[2] = vhat_jp1[2] - zhat_jp1[2];
      for (nap = 0; nap <= 994; nap += 2) {
        r = _mm_loadu_pd(&zhat_jp1[nap + 3]);
        _mm_storeu_pd(&c[nap + 3], _mm_mul_pd(r, _mm_set1_pd(-1.0)));
      }
      c[999] = -zhat_jp1[999];
      for (i = 0; i <= 998; i += 2) {
        r = _mm_loadu_pd(&zhat_jp1[i]);
        r1 = _mm_loadu_pd(&zhat_j[i]);
        _mm_storeu_pd(&zhat_j[i], _mm_sub_pd(r, r1));
      }
      n_t = (ptrdiff_t)1000;
      incx_t = (ptrdiff_t)1;
      incy_t = (ptrdiff_t)1;
      scale = ddot(&n_t, &c[0], &incx_t, &zhat_j[0], &incy_t);
      if (((scale <= 0.0) && (j_restart == -1.0)) ||
          (b_mod((real_T)j + 1.0, j_restart) != 0.0)) {
        j_r++;
      } else {
        j_r = 1.0;
        /* fprintf("Restarting at j = %d\n", j) */
      }
      /*  zhat_inf_dj2 = zhat_inf_djp1; */
      /*  what_inf_dj2 = what_inf_djp1; */
      /*     %% Update values for next iteration */
      memcpy(&zhat_j[0], &zhat_jp1[0], 1000U * sizeof(real_T));
      what_j[0] = what_star[0];
      vhat_j[0] = vhat_jp1[0];
      what_j[1] = what_star[1];
      vhat_j[1] = vhat_jp1[1];
      what_j[2] = what_star[2];
      vhat_j[2] = vhat_jp1[2];
      j++;
    }
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  emxFree_int32_T(sp, &a_rowidx);
  emxFree_int32_T(sp, &a_colidx);
  emxFree_real_T(sp, &a_d);
  /*  Recover original primal variables */
  st.site = &g_emlrtRSI;
  for (i = 0; i <= 999998; i += 2) {
    _mm_storeu_pd(&SD->f0.S_z[i],
                  _mm_mul_pd(_mm_set1_pd(S_z), _mm_loadu_pd(&L_inv[i])));
  }
  b_st.site = &n_emlrtRSI;
  mtimes(SD->f0.S_z, zhat_jp1, z_star);
  for (i = 0; i <= 998; i += 2) {
    r = _mm_loadu_pd(&z_star[i]);
    _mm_storeu_pd(&z_star[i], _mm_add_pd(r, _mm_set1_pd(c_z)));
  }
  /*  Retain transformed dual variable */
  st.site = &h_emlrtRSI;
  scale = toc(&st, expl_temp.tv_sec, expl_temp.tv_nsec);
  st.site = &i_emlrtRSI;
  b_st.site = &ob_emlrtRSI;
  y = NULL;
  propValues = emlrtCreateCharArray(2, &iv[0]);
  emlrtInitCharArrayR2013a(&b_st, 7, propValues, &u[0]);
  emlrtAssign(&y, propValues);
  b_y = NULL;
  propValues = emlrtCreateDoubleScalar(1.0);
  emlrtAssign(&b_y, propValues);
  c_y = NULL;
  propValues = emlrtCreateCharArray(2, &iv1[0]);
  emlrtInitCharArrayR2013a(&b_st, 20, propValues, &b_u[0]);
  emlrtAssign(&c_y, propValues);
  d_y = NULL;
  propValues = emlrtCreateDoubleScalar(scale * 1000.0);
  emlrtAssign(&d_y, propValues);
  c_st.site = &qb_emlrtRSI;
  scale = emlrt_marshallIn(&c_st, feval(&c_st, y, b_y, c_y, d_y, &b_emlrtMCI),
                           "<output of feval>");
  e_y = NULL;
  propValues = emlrtCreateDoubleScalar(scale);
  emlrtAssign(&e_y, propValues);
  emlrtDisplayR2012b(e_y, "ans", &emlrtRTEI, (emlrtCTX)sp);
  /*   */
  /*  figure */
  /*  plot(aff_violation) */
  /*  yscale("log") */
  /*   */
  /*  figure */
  /*  zhis =  zhis; */
  /*  scatter(zhis(1, :), zhis(2, :)) */
  /*  axis equal */
  /*  grid on */
  /*  %  */
  /*  % figure */
  /*  % zhis = zhis; */
  /*  % plot(vecnorm(zhis - zhis(:, end))) */
  /*  % yscale("log") */
  /*  % grid on */
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (PIPG_orig_precond_Dexplicit.c) */
