/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_custom_power_iteration_noops_mex.c
 *
 * Code generation for function '_coder_custom_power_iteration_noops_mex'
 *
 */

/* Include files */
#include "_coder_custom_power_iteration_noops_mex.h"
#include "_coder_custom_power_iteration_noops_api.h"
#include "custom_power_iteration_noops_data.h"
#include "custom_power_iteration_noops_initialize.h"
#include "custom_power_iteration_noops_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void custom_power_iteration_noops_mexFunction(int32_T nlhs, mxArray *plhs[1],
                                              int32_T nrhs,
                                              const mxArray *prhs[11])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 11) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 11, 4,
                        28, "custom_power_iteration_noops");
  }
  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 28,
                        "custom_power_iteration_noops");
  }
  /* Call the function. */
  e_custom_power_iteration_noops_(prhs, &outputs);
  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, &plhs[0], &outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&custom_power_iteration_noops_atexit);
  custom_power_iteration_noops_initialize();
  custom_power_iteration_noops_mexFunction(nlhs, plhs, nrhs, prhs);
  custom_power_iteration_noops_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_custom_power_iteration_noops_mex.c) */
