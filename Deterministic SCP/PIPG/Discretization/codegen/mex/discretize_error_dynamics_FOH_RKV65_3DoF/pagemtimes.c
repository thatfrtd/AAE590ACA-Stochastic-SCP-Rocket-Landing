/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * pagemtimes.c
 *
 * Code generation for function 'pagemtimes'
 *
 */

/* Include files */
#include "pagemtimes.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_data.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_emxutil.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo u_emlrtRSI = {
    25,           /* lineNo */
    "pagemtimes", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\pagemtimes.m" /* pathName
                                                                           */
};

static emlrtRSInfo v_emlrtRSI = {
    61,           /* lineNo */
    "pagemtimes", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\pagemtimes.m" /* pathName
                                                                           */
};

static emlrtRSInfo w_emlrtRSI = {
    66,           /* lineNo */
    "pagemtimes", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\pagemtimes.m" /* pathName
                                                                           */
};

static emlrtRSInfo x_emlrtRSI = {
    67,           /* lineNo */
    "pagemtimes", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\pagemtimes.m" /* pathName
                                                                           */
};

static emlrtRSInfo ab_emlrtRSI = {
    129,                          /* lineNo */
    "computeMultiplicationIndex", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\pagemtimes.m" /* pathName
                                                                           */
};

static emlrtRSInfo bb_emlrtRSI = {
    33,        /* lineNo */
    "idivide", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\idivide.m" /* pathName
                                                                        */
};

static emlrtRSInfo cb_emlrtRSI = {
    81,            /* lineNo */
    "idivide_fix", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\idivide.m" /* pathName
                                                                        */
};

static emlrtRTEInfo g_emlrtRTEI = {
    132,           /* lineNo */
    9,             /* colNo */
    "eml_idivide", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\idivide.m" /* pName
                                                                        */
};

static emlrtRTEInfo h_emlrtRTEI = {
    188,                 /* lineNo */
    13,                  /* colNo */
    "processArraySizes", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\pagemtimes.m" /* pName
                                                                           */
};

static emlrtRTEInfo fd_emlrtRTEI = {
    56,           /* lineNo */
    24,           /* colNo */
    "pagemtimes", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\pagemtimes.m" /* pName
                                                                           */
};

static emlrtRSInfo sb_emlrtRSI = {
    103,      /* lineNo */
    "intmod", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\eml\\+coder\\+internal\\+scalar\\mod."
    "m" /* pathName */
};

/* Function Declarations */
static int32_T div_s32(const emlrtStack *sp, int32_T numerator,
                       int32_T denominator);

/* Function Definitions */
static int32_T div_s32(const emlrtStack *sp, int32_T numerator,
                       int32_T denominator)
{
  int32_T quotient;
  if (denominator == 0) {
    emlrtDivisionByZeroErrorR2012b(NULL, (emlrtConstCTX)sp);
  } else {
    uint32_T u;
    uint32_T u1;
    if (numerator < 0) {
      u = ~(uint32_T)numerator + 1U;
    } else {
      u = (uint32_T)numerator;
    }
    if (denominator < 0) {
      u1 = ~(uint32_T)denominator + 1U;
    } else {
      u1 = (uint32_T)denominator;
    }
    u /= u1;
    if ((numerator < 0) != (denominator < 0)) {
      quotient = -(int32_T)u;
    } else {
      quotient = (int32_T)u;
    }
  }
  return quotient;
}

void b_pagemtimes(const emlrtStack *sp, const emxArray_real_T *varargin_1,
                  const emxArray_real_T *varargin_2, emxArray_real_T *Z)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  const real_T *varargin_1_data;
  const real_T *varargin_2_data;
  real_T *Z_data;
  int32_T b_i;
  int32_T i;
  int32_T i1;
  int32_T ii;
  int32_T j;
  int32_T k;
  int32_T outputSize_idx_2;
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  varargin_2_data = varargin_2->data;
  varargin_1_data = varargin_1->data;
  st.site = &u_emlrtRSI;
  if (varargin_1->size[2] != 1) {
    if ((varargin_2->size[2] != 1) &&
        (varargin_1->size[2] != varargin_2->size[2])) {
      emlrtErrorWithMessageIdR2018a(
          &st, &h_emlrtRTEI, "MATLAB:pagemtimes:mismatchImplicitExpansion",
          "MATLAB:pagemtimes:mismatchImplicitExpansion", 0);
    }
    outputSize_idx_2 = varargin_1->size[2];
  } else {
    outputSize_idx_2 = varargin_2->size[2];
  }
  i = Z->size[0] * Z->size[1] * Z->size[2];
  Z->size[0] = 7;
  Z->size[1] = 2;
  Z->size[2] = outputSize_idx_2;
  emxEnsureCapacity_real_T(sp, Z, i, &fd_emlrtRTEI);
  Z_data = Z->data;
  p = (varargin_1->size[2] == varargin_2->size[2]);
  st.site = &v_emlrtRSI;
  if (outputSize_idx_2 > 2147483646) {
    b_st.site = &y_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }
  for (ii = 0; ii < outputSize_idx_2; ii++) {
    real_T C[14];
    int32_T X_idx;
    int32_T Y_idx;
    if (p) {
      X_idx = ii;
      Y_idx = ii;
    } else {
      st.site = &w_emlrtRSI;
      X_idx = 0;
      if (varargin_1->size[2] != 1) {
        if (outputSize_idx_2 == 0) {
          X_idx = ii;
        } else {
          b_st.site = &sb_emlrtRSI;
          X_idx = ii - div_s32(&b_st, ii, outputSize_idx_2) * outputSize_idx_2;
        }
      }
      b_st.site = &ab_emlrtRSI;
      c_st.site = &bb_emlrtRSI;
      d_st.site = &cb_emlrtRSI;
      if (outputSize_idx_2 == 0) {
        emlrtErrorWithMessageIdR2018a(&d_st, &g_emlrtRTEI,
                                      "Coder:toolbox:idivide_divideByZero",
                                      "Coder:toolbox:idivide_divideByZero", 0);
      }
      st.site = &x_emlrtRSI;
      Y_idx = 0;
      if (varargin_2->size[2] != 1) {
        b_st.site = &sb_emlrtRSI;
        Y_idx = ii - div_s32(&b_st, ii, outputSize_idx_2) * outputSize_idx_2;
      }
      b_st.site = &ab_emlrtRSI;
      c_st.site = &bb_emlrtRSI;
    }
    for (j = 0; j < 2; j++) {
      int32_T coffset_tmp;
      coffset_tmp = j * 7;
      for (b_i = 0; b_i < 7; b_i++) {
        real_T s;
        s = 0.0;
        for (k = 0; k < 7; k++) {
          i = k * 7 + b_i;
          i1 = coffset_tmp + k;
          s += varargin_1_data[(i % 7 + 7 * (i / 7)) + 49 * X_idx] *
               varargin_2_data[(i1 % 7 + 7 * (i1 / 7)) + 14 * Y_idx];
        }
        C[coffset_tmp + b_i] = s;
      }
    }
    for (i = 0; i < 2; i++) {
      for (i1 = 0; i1 < 7; i1++) {
        Z_data[(i1 + Z->size[0] * i) + Z->size[0] * Z->size[1] * ii] =
            C[i1 + 7 * i];
      }
    }
  }
}

void c_pagemtimes(const emlrtStack *sp, const emxArray_real_T *varargin_1,
                  const emxArray_real_T *varargin_2, emxArray_real_T *Z)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  const real_T *varargin_1_data;
  const real_T *varargin_2_data;
  real_T *Z_data;
  int32_T b_i;
  int32_T i;
  int32_T ii;
  int32_T k;
  int32_T outputSize_idx_2;
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  varargin_2_data = varargin_2->data;
  varargin_1_data = varargin_1->data;
  st.site = &u_emlrtRSI;
  if (varargin_1->size[2] != 1) {
    if ((varargin_2->size[2] != 1) &&
        (varargin_1->size[2] != varargin_2->size[2])) {
      emlrtErrorWithMessageIdR2018a(
          &st, &h_emlrtRTEI, "MATLAB:pagemtimes:mismatchImplicitExpansion",
          "MATLAB:pagemtimes:mismatchImplicitExpansion", 0);
    }
    outputSize_idx_2 = varargin_1->size[2];
  } else {
    outputSize_idx_2 = varargin_2->size[2];
  }
  i = Z->size[0] * Z->size[1] * Z->size[2];
  Z->size[0] = 7;
  Z->size[1] = 1;
  Z->size[2] = outputSize_idx_2;
  emxEnsureCapacity_real_T(sp, Z, i, &fd_emlrtRTEI);
  Z_data = Z->data;
  p = (varargin_1->size[2] == varargin_2->size[2]);
  st.site = &v_emlrtRSI;
  if (outputSize_idx_2 > 2147483646) {
    b_st.site = &y_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }
  for (ii = 0; ii < outputSize_idx_2; ii++) {
    real_T C[7];
    int32_T X_idx;
    int32_T Y_idx;
    if (p) {
      X_idx = ii;
      Y_idx = ii;
    } else {
      st.site = &w_emlrtRSI;
      X_idx = 0;
      if (varargin_1->size[2] != 1) {
        if (outputSize_idx_2 == 0) {
          X_idx = ii;
        } else {
          b_st.site = &sb_emlrtRSI;
          X_idx = ii - div_s32(&b_st, ii, outputSize_idx_2) * outputSize_idx_2;
        }
      }
      b_st.site = &ab_emlrtRSI;
      c_st.site = &bb_emlrtRSI;
      d_st.site = &cb_emlrtRSI;
      if (outputSize_idx_2 == 0) {
        emlrtErrorWithMessageIdR2018a(&d_st, &g_emlrtRTEI,
                                      "Coder:toolbox:idivide_divideByZero",
                                      "Coder:toolbox:idivide_divideByZero", 0);
      }
      st.site = &x_emlrtRSI;
      Y_idx = 0;
      if (varargin_2->size[2] != 1) {
        b_st.site = &sb_emlrtRSI;
        Y_idx = ii - div_s32(&b_st, ii, outputSize_idx_2) * outputSize_idx_2;
      }
      b_st.site = &ab_emlrtRSI;
      c_st.site = &bb_emlrtRSI;
    }
    for (b_i = 0; b_i < 7; b_i++) {
      real_T s;
      s = 0.0;
      for (k = 0; k < 7; k++) {
        i = k * 7 + b_i;
        s += varargin_1_data[(i % 7 + 7 * (i / 7)) + 49 * X_idx] *
             varargin_2_data[k + 7 * Y_idx];
      }
      C[b_i] = s;
    }
    for (i = 0; i < 7; i++) {
      Z_data[i + Z->size[0] * Z->size[1] * ii] = C[i];
    }
  }
}

void d_pagemtimes(const emlrtStack *sp, const emxArray_real_T *varargin_1,
                  const real_T varargin_2_data[],
                  const int32_T varargin_2_size[3], emxArray_real_T *Z)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  const real_T *varargin_1_data;
  real_T *Z_data;
  int32_T X_idx;
  int32_T i;
  int32_T ii;
  int32_T outputSize_idx_2;
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  varargin_1_data = varargin_1->data;
  st.site = &u_emlrtRSI;
  if (varargin_1->size[2] != 1) {
    if ((varargin_2_size[2] != 1) &&
        (varargin_1->size[2] != varargin_2_size[2])) {
      emlrtErrorWithMessageIdR2018a(
          &st, &h_emlrtRTEI, "MATLAB:pagemtimes:mismatchImplicitExpansion",
          "MATLAB:pagemtimes:mismatchImplicitExpansion", 0);
    }
    outputSize_idx_2 = varargin_1->size[2];
  } else {
    outputSize_idx_2 = varargin_2_size[2];
  }
  X_idx = Z->size[0] * Z->size[1] * Z->size[2];
  Z->size[0] = 7;
  Z->size[1] = 1;
  Z->size[2] = outputSize_idx_2;
  emxEnsureCapacity_real_T(sp, Z, X_idx, &fd_emlrtRTEI);
  Z_data = Z->data;
  p = (varargin_1->size[2] == varargin_2_size[2]);
  st.site = &v_emlrtRSI;
  if (outputSize_idx_2 > 2147483646) {
    b_st.site = &y_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }
  for (ii = 0; ii < outputSize_idx_2; ii++) {
    real_T C[7];
    int32_T Y_idx;
    if (p) {
      X_idx = ii;
      Y_idx = ii;
    } else {
      st.site = &w_emlrtRSI;
      X_idx = 0;
      if (varargin_1->size[2] != 1) {
        if (outputSize_idx_2 == 0) {
          X_idx = ii;
        } else {
          b_st.site = &sb_emlrtRSI;
          X_idx = ii - div_s32(&b_st, ii, outputSize_idx_2) * outputSize_idx_2;
        }
      }
      b_st.site = &ab_emlrtRSI;
      c_st.site = &bb_emlrtRSI;
      d_st.site = &cb_emlrtRSI;
      if (outputSize_idx_2 == 0) {
        emlrtErrorWithMessageIdR2018a(&d_st, &g_emlrtRTEI,
                                      "Coder:toolbox:idivide_divideByZero",
                                      "Coder:toolbox:idivide_divideByZero", 0);
      }
      st.site = &x_emlrtRSI;
      Y_idx = 0;
      if (varargin_2_size[2] != 1) {
        b_st.site = &sb_emlrtRSI;
        Y_idx = ii - div_s32(&b_st, ii, outputSize_idx_2) * outputSize_idx_2;
      }
      b_st.site = &ab_emlrtRSI;
      c_st.site = &bb_emlrtRSI;
    }
    for (i = 0; i < 7; i++) {
      C[i] = varargin_1_data[i % 7 + 14 * X_idx] * varargin_2_data[2 * Y_idx] +
             varargin_1_data[((i + 7) % 7 + 14 * X_idx) + 7] *
                 varargin_2_data[2 * Y_idx + 1];
    }
    for (X_idx = 0; X_idx < 7; X_idx++) {
      Z_data[X_idx + Z->size[0] * Z->size[1] * ii] = C[X_idx];
    }
  }
}

void pagemtimes(const emlrtStack *sp, const emxArray_real_T *varargin_1,
                const emxArray_real_T *varargin_2, emxArray_real_T *Z)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  const real_T *varargin_1_data;
  const real_T *varargin_2_data;
  real_T *Z_data;
  int32_T b_i;
  int32_T i;
  int32_T ii;
  int32_T j;
  int32_T k;
  int32_T outputSize_idx_2;
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  varargin_2_data = varargin_2->data;
  varargin_1_data = varargin_1->data;
  st.site = &u_emlrtRSI;
  if (varargin_1->size[2] != 1) {
    if ((varargin_2->size[2] != 1) &&
        (varargin_1->size[2] != varargin_2->size[2])) {
      emlrtErrorWithMessageIdR2018a(
          &st, &h_emlrtRTEI, "MATLAB:pagemtimes:mismatchImplicitExpansion",
          "MATLAB:pagemtimes:mismatchImplicitExpansion", 0);
    }
    outputSize_idx_2 = varargin_1->size[2];
  } else {
    outputSize_idx_2 = varargin_2->size[2];
  }
  i = Z->size[0] * Z->size[1] * Z->size[2];
  Z->size[0] = 7;
  Z->size[1] = 7;
  Z->size[2] = outputSize_idx_2;
  emxEnsureCapacity_real_T(sp, Z, i, &fd_emlrtRTEI);
  Z_data = Z->data;
  p = (varargin_1->size[2] == varargin_2->size[2]);
  st.site = &v_emlrtRSI;
  if (outputSize_idx_2 > 2147483646) {
    b_st.site = &y_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }
  for (ii = 0; ii < outputSize_idx_2; ii++) {
    real_T C[49];
    int32_T X_idx;
    int32_T Y_idx;
    if (p) {
      X_idx = ii;
      Y_idx = ii;
    } else {
      st.site = &w_emlrtRSI;
      X_idx = 0;
      if (varargin_1->size[2] != 1) {
        if (outputSize_idx_2 == 0) {
          X_idx = ii;
        } else {
          b_st.site = &sb_emlrtRSI;
          X_idx = ii - div_s32(&b_st, ii, outputSize_idx_2) * outputSize_idx_2;
        }
      }
      b_st.site = &ab_emlrtRSI;
      c_st.site = &bb_emlrtRSI;
      d_st.site = &cb_emlrtRSI;
      if (outputSize_idx_2 == 0) {
        emlrtErrorWithMessageIdR2018a(&d_st, &g_emlrtRTEI,
                                      "Coder:toolbox:idivide_divideByZero",
                                      "Coder:toolbox:idivide_divideByZero", 0);
      }
      st.site = &x_emlrtRSI;
      Y_idx = 0;
      if (varargin_2->size[2] != 1) {
        b_st.site = &sb_emlrtRSI;
        Y_idx = ii - div_s32(&b_st, ii, outputSize_idx_2) * outputSize_idx_2;
      }
      b_st.site = &ab_emlrtRSI;
      c_st.site = &bb_emlrtRSI;
    }
    for (j = 0; j < 7; j++) {
      int32_T coffset_tmp;
      coffset_tmp = j * 7;
      for (b_i = 0; b_i < 7; b_i++) {
        real_T s;
        s = 0.0;
        for (k = 0; k < 7; k++) {
          int32_T i1;
          i = k * 7 + b_i;
          i1 = coffset_tmp + k;
          s += varargin_1_data[(i % 7 + 7 * (i / 7)) + 49 * X_idx] *
               varargin_2_data[(i1 % 7 + 7 * (i1 / 7)) + 49 * Y_idx];
        }
        C[coffset_tmp + b_i] = s;
      }
    }
    for (i = 0; i < 49; i++) {
      Z_data[i + ii * 49] = C[i];
    }
  }
}

/* End of code generation (pagemtimes.c) */
