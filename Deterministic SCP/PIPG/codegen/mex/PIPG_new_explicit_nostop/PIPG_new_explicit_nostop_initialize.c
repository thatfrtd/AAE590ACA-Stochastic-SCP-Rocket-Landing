/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_new_explicit_nostop_initialize.c
 *
 * Code generation for function 'PIPG_new_explicit_nostop_initialize'
 *
 */

/* Include files */
#include "PIPG_new_explicit_nostop_initialize.h"
#include "PIPG_new_explicit_nostop_data.h"
#include "_coder_PIPG_new_explicit_nostop_mex.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void PIPG_new_explicit_nostop_once(void);

/* Function Definitions */
static void PIPG_new_explicit_nostop_once(void)
{
  static const int32_T lineInfo[38] = {9,  10, 11, 12, 13, 14, 15,  16, 18, 20,
                                       23, 26, 29, 30, 33, 34, 37,  39, 40, 42,
                                       44, 46, 47, 50, 53, 56, 59,  62, 64, 66,
                                       68, 69, 72, 75, 77, 78, 101, 103};
  static const int32_T b_lineInfo[6] = {108, 110, 111, 113, 115, 116};
  static const int32_T c_lineInfo[3] = {122, 124, 125};
  mex_InitInfAndNan();
  singleton_project_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\PIPG_new_explicit_nostop.m>singleton_project(codegen)";
  ball_project_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\PIPG_new_explicit_nostop.m>ball_project(codegen)";
  c_PIPG_new_explicit_nostop_comp =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\PIPG_new_explicit_nostop.m>PIPG_new_explicit_nostop(codegen)";
  isMexOutdated = emlrtProfilerCheckMEXOutdated();
  emlrtProfilerRegisterMEXFcn(
      (char_T *)c_PIPG_new_explicit_nostop_comp,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-"
                "Landing\\Deterministic SCP\\P"
                "IPG\\PIPG_new_explicit_nostop.m",
      (char_T *)"PIPG_new_explicit_nostop", 38, &lineInfo[0], isMexOutdated);
  emlrtProfilerRegisterMEXFcn(
      (char_T *)ball_project_complete_name,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-"
                "Landing\\Deterministic SCP\\P"
                "IPG\\PIPG_new_explicit_nostop.m",
      (char_T *)"PIPG_new_explicit_nostop>ball_project", 6, &b_lineInfo[0],
      isMexOutdated);
  emlrtProfilerRegisterMEXFcn(
      (char_T *)singleton_project_complete_name,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-"
                "Landing\\Deterministic SCP\\P"
                "IPG\\PIPG_new_explicit_nostop.m",
      (char_T *)"PIPG_new_explicit_nostop>singleton_project", 3, &c_lineInfo[0],
      isMexOutdated);
}

void PIPG_new_explicit_nostop_initialize(void)
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
    PIPG_new_explicit_nostop_once();
  }
  emlrtCheckProfilerStatus();
}

/* End of code generation (PIPG_new_explicit_nostop_initialize.c) */
