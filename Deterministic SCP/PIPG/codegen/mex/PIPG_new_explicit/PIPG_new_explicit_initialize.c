/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_new_explicit_initialize.c
 *
 * Code generation for function 'PIPG_new_explicit_initialize'
 *
 */

/* Include files */
#include "PIPG_new_explicit_initialize.h"
#include "PIPG_new_explicit_data.h"
#include "_coder_PIPG_new_explicit_mex.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void PIPG_new_explicit_once(void);

/* Function Definitions */
static void PIPG_new_explicit_once(void)
{
  static const int32_T lineInfo[67] = {
      9,  10,  11,  12,  13,  14,  15,  16,  17,  19,  21, 24, 27, 30,
      31, 34,  35,  38,  40,  41,  43,  45,  47,  48,  51, 54, 57, 60,
      63, 64,  65,  66,  68,  69,  70,  72,  73,  74,  75, 76, 77, 78,
      79, 81,  82,  83,  84,  85,  86,  88,  89,  90,  91, 92, 93, 96,
      98, 100, 102, 103, 106, 109, 111, 112, 114, 135, 137};
  static const int32_T d_lineInfo[12] = {5,  6,  7,  9,  10, 11,
                                         13, 14, 15, 16, 17, 19};
  static const int32_T b_lineInfo[10] = {23, 25, 26, 27, 29,
                                         31, 32, 34, 36, 37};
  static const int32_T c_lineInfo[9] = {24, 25, 26, 28, 30, 31, 33, 35, 36};
  mex_InitInfAndNan();
  stopping_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\stopping.m>stopping(codegen)";
  c_Singletonz_project_complete_n =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\Convex Cones\\Singletonz.m>Singletonz.project(codegen)";
  Ballz_project_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\Convex Cones\\Ballz.m>Ballz.project(codegen)";
  PIPG_new_explicit_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\PIPG_new_explicit.m>PIPG_new_explicit(codegen)";
  isMexOutdated = emlrtProfilerCheckMEXOutdated();
  emlrtProfilerRegisterMEXFcn(
      (char_T *)PIPG_new_explicit_complete_name,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-"
                "Landing\\Deterministic SCP\\P"
                "IPG\\PIPG_new_explicit.m",
      (char_T *)"PIPG_new_explicit", 67, &lineInfo[0], isMexOutdated);
  emlrtProfilerRegisterMEXFcn(
      (char_T *)Ballz_project_complete_name,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-"
                "Landing\\Deterministic SCP\\P"
                "IPG\\Convex Cones\\Ballz.m",
      (char_T *)"Ballz>Ballz.project", 10, &b_lineInfo[0], isMexOutdated);
  emlrtProfilerRegisterMEXFcn(
      (char_T *)c_Singletonz_project_complete_n,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-"
                "Landing\\Deterministic SCP\\P"
                "IPG\\Convex Cones\\Singletonz.m",
      (char_T *)"Singletonz>Singletonz.project", 9, &c_lineInfo[0],
      isMexOutdated);
  emlrtProfilerRegisterMEXFcn(
      (char_T *)stopping_complete_name,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-"
                "Landing\\Deterministic SCP\\P"
                "IPG\\stopping.m",
      (char_T *)"stopping", 12, &d_lineInfo[0], isMexOutdated);
}

void PIPG_new_explicit_initialize(void)
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
    PIPG_new_explicit_once();
  }
  emlrtCheckProfilerStatus();
}

/* End of code generation (PIPG_new_explicit_initialize.c) */
