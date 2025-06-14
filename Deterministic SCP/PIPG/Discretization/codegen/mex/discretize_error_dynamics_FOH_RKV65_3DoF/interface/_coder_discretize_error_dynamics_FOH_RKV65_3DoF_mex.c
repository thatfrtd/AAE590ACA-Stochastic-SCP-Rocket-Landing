/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_discretize_error_dynamics_FOH_RKV65_3DoF_mex.c
 *
 * Code generation for function
 * '_coder_discretize_error_dynamics_FOH_RKV65_3DoF_mex'
 *
 */

/* Include files */
#include "_coder_discretize_error_dynamics_FOH_RKV65_3DoF_mex.h"
#include "_coder_discretize_error_dynamics_FOH_RKV65_3DoF_api.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_data.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_initialize.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void discretize_error_dynamics_FOH_RKV65_3DoF_mexFunction(
    int32_T nlhs, mxArray *plhs[6], int32_T nrhs, const mxArray *prhs[9])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[6];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 9) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 9, 4,
                        40, "discretize_error_dynamics_FOH_RKV65_3DoF");
  }
  if (nlhs > 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 40,
                        "discretize_error_dynamics_FOH_RKV65_3DoF");
  }
  /* Call the function. */
  e_discretize_error_dynamics_FOH(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&discretize_error_dynamics_FOH_RKV65_3DoF_atexit);
  discretize_error_dynamics_FOH_RKV65_3DoF_initialize();
  discretize_error_dynamics_FOH_RKV65_3DoF_mexFunction(nlhs, plhs, nrhs, prhs);
  discretize_error_dynamics_FOH_RKV65_3DoF_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation
 * (_coder_discretize_error_dynamics_FOH_RKV65_3DoF_mex.c) */
