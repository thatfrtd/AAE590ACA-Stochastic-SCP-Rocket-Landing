/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PIPG_new_explicit_mex.c
 *
 * Code generation for function '_coder_PIPG_new_explicit_mex'
 *
 */

/* Include files */
#include "_coder_PIPG_new_explicit_mex.h"
#include "PIPG_new_explicit_data.h"
#include "PIPG_new_explicit_initialize.h"
#include "PIPG_new_explicit_terminate.h"
#include "_coder_PIPG_new_explicit_api.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void PIPG_new_explicit_mexFunction(int32_T nlhs, mxArray *plhs[3], int32_T nrhs,
                                   const mxArray *prhs[18])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[3];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 18) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 18, 4,
                        17, "PIPG_new_explicit");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 17,
                        "PIPG_new_explicit");
  }
  /* Call the function. */
  PIPG_new_explicit_api(prhs, nlhs, outputs);
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
  mexAtExit(&PIPG_new_explicit_atexit);
  PIPG_new_explicit_initialize();
  PIPG_new_explicit_mexFunction(nlhs, plhs, nrhs, prhs);
  PIPG_new_explicit_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_PIPG_new_explicit_mex.c) */
