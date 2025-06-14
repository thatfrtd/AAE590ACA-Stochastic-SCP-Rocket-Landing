/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * custom_power_iteration_noops_initialize.c
 *
 * Code generation for function 'custom_power_iteration_noops_initialize'
 *
 */

/* Include files */
#include "custom_power_iteration_noops_initialize.h"
#include "_coder_custom_power_iteration_noops_mex.h"
#include "custom_power_iteration_noops_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void c_custom_power_iteration_noops_(void);

/* Function Definitions */
static void c_custom_power_iteration_noops_(void)
{
  static const int32_T lineInfo[48] = {
      5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21,
      22, 24, 26, 27, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40, 42,
      43, 44, 45, 47, 48, 49, 51, 52, 54, 55, 56, 57, 58, 59, 61, 63};
  mex_InitInfAndNan();
  d_custom_power_iteration_noops_ =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\custom_power_iteration_noops.m>custom_power_iteration_noops("
      "codegen)";
  isMexOutdated = emlrtProfilerCheckMEXOutdated();
  emlrtProfilerRegisterMEXFcn(
      (char_T *)d_custom_power_iteration_noops_,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-"
                "Landing\\Deterministic SCP\\P"
                "IPG\\custom_power_iteration_noops.m",
      (char_T *)"custom_power_iteration_noops", 48, &lineInfo[0],
      isMexOutdated);
}

void custom_power_iteration_noops_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    c_custom_power_iteration_noops_();
  }
  emlrtCheckProfilerStatus();
}

/* End of code generation (custom_power_iteration_noops_initialize.c) */
