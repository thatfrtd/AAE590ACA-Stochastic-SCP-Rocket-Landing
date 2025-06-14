/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * zero_if_empty.c
 *
 * Code generation for function 'zero_if_empty'
 *
 */

/* Include files */
#include "zero_if_empty.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_data.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_emxutil.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo gd_emlrtRTEI = {
    6,               /* lineNo */
    1,               /* colNo */
    "zero_if_empty", /* fName */
    "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
    "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Helper Functions\\Ge"
    "neral\\zero_if_empty.m" /* pName */
};

/* Function Definitions */
void zero_if_empty(const emlrtStack *sp, const emxArray_real_T *v,
                   emxArray_real_T *b_v)
{
  real_T *v_data;
  int32_T i;
  int32_T loop_ub_tmp;
  uint32_T new_size_idx_2;
  emlrtMEXProfilingFunctionEntry((char_T *)zero_if_empty_complete_name,
                                 isMexOutdated);
  /* ZERO_IF_EMPTY Summary of this function goes here */
  /*    Detailed explanation goes here */
  emlrtMEXProfilingStatement(1, isMexOutdated);
  new_size_idx_2 = (uint32_T)v->size[2];
  emlrtMEXProfilingStatement(2, isMexOutdated);
  if (v->size[2] == 0) {
    new_size_idx_2 = 1U;
  }
  emlrtMEXProfilingStatement(3, isMexOutdated);
  i = b_v->size[0] * b_v->size[1] * b_v->size[2];
  b_v->size[0] = 7;
  b_v->size[1] = 1;
  b_v->size[2] = (int32_T)new_size_idx_2;
  emxEnsureCapacity_real_T(sp, b_v, i, &gd_emlrtRTEI);
  v_data = b_v->data;
  loop_ub_tmp = 7 * (int32_T)new_size_idx_2;
  for (i = 0; i < loop_ub_tmp; i++) {
    v_data[i] = 0.0;
  }
  emlrtMEXProfilingStatement(4, isMexOutdated);
  emlrtMEXProfilingFunctionExit(isMexOutdated);
}

/* End of code generation (zero_if_empty.c) */
