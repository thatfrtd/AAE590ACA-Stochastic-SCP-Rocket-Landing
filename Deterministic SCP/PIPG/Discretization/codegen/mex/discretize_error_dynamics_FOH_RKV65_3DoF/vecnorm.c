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
#include "discretize_error_dynamics_FOH_RKV65_3DoF_data.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_emxutil.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo eb_emlrtRSI = {
    51,        /* lineNo */
    "vecnorm", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\matfun\\vecnorm.m" /* pathName
                                                                           */
};

static emlrtRSInfo fb_emlrtRSI = {
    52,        /* lineNo */
    "vecnorm", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\matfun\\vecnorm.m" /* pathName
                                                                           */
};

static emlrtRSInfo gb_emlrtRSI = {
    23,      /* lineNo */
    "xnrm2", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xnrm2."
    "m" /* pathName */
};

static emlrtRSInfo hb_emlrtRSI = {
    38,      /* lineNo */
    "xnrm2", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\eml\\+coder\\+internal\\+"
    "refblas\\xnrm2.m" /* pathName */
};

static emlrtRTEInfo hd_emlrtRTEI = {
    50,        /* lineNo */
    24,        /* colNo */
    "vecnorm", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\matfun\\vecnorm.m" /* pName
                                                                           */
};

/* Function Definitions */
void vecnorm(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  const real_T *x_data;
  real_T *y_data;
  int32_T ix0;
  int32_T j;
  int32_T k;
  int32_T ncols_tmp;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  x_data = x->data;
  ncols_tmp = x->size[1];
  ix0 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity_real_T(sp, y, ix0, &hd_emlrtRTEI);
  y_data = y->data;
  st.site = &eb_emlrtRSI;
  if (x->size[1] > 2147483646) {
    b_st.site = &y_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }
  for (j = 0; j < ncols_tmp; j++) {
    real_T b_y;
    real_T scale;
    st.site = &fb_emlrtRSI;
    ix0 = j * 7 + 7;
    b_st.site = &gb_emlrtRSI;
    b_y = 0.0;
    scale = 3.3121686421112381E-170;
    c_st.site = &hb_emlrtRSI;
    if ((ix0 - 6 <= ix0) && (ix0 > 2147483646)) {
      d_st.site = &y_emlrtRSI;
      check_forloop_overflow_error(&d_st);
    }
    for (k = ix0 - 6; k <= ix0; k++) {
      real_T absxk;
      absxk = muDoubleScalarAbs(x_data[k - 1]);
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
    y_data[j] = scale * muDoubleScalarSqrt(b_y);
  }
}

/* End of code generation (vecnorm.c) */
