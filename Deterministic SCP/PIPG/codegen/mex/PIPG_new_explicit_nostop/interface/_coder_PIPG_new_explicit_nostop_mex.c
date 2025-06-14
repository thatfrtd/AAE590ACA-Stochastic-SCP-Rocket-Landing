/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PIPG_new_explicit_nostop_mex.c
 *
 * Code generation for function '_coder_PIPG_new_explicit_nostop_mex'
 *
 */

/* Include files */
#include "_coder_PIPG_new_explicit_nostop_mex.h"
#include "PIPG_new_explicit_nostop_data.h"
#include "PIPG_new_explicit_nostop_initialize.h"
#include "PIPG_new_explicit_nostop_terminate.h"
#include "_coder_PIPG_new_explicit_nostop_api.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void PIPG_new_explicit_nostop_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                          int32_T nrhs, const mxArray *prhs[16])
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
  if (nrhs != 16) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 16, 4,
                        24, "PIPG_new_explicit_nostop");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 24,
                        "PIPG_new_explicit_nostop");
  }
  /* Call the function. */
  PIPG_new_explicit_nostop_api(prhs, nlhs, outputs);
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
  mexAtExit(&PIPG_new_explicit_nostop_atexit);
  PIPG_new_explicit_nostop_initialize();
  PIPG_new_explicit_nostop_mexFunction(nlhs, plhs, nrhs, prhs);
  PIPG_new_explicit_nostop_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_PIPG_new_explicit_nostop_mex.c) */
