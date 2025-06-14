/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PIPG_orig_precond_Dexplicit_mex.c
 *
 * Code generation for function '_coder_PIPG_orig_precond_Dexplicit_mex'
 *
 */

/* Include files */
#include "_coder_PIPG_orig_precond_Dexplicit_mex.h"
#include "PIPG_orig_precond_Dexplicit_data.h"
#include "PIPG_orig_precond_Dexplicit_initialize.h"
#include "PIPG_orig_precond_Dexplicit_terminate.h"
#include "PIPG_orig_precond_Dexplicit_types.h"
#include "_coder_PIPG_orig_precond_Dexplicit_api.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void PIPG_orig_precond_Dexplicit_mexFunction(
    c_PIPG_orig_precond_DexplicitSt *SD, int32_T nlhs, mxArray *plhs[3],
    int32_T nrhs, const mxArray *prhs[20])
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
  if (nrhs != 20) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 20, 4,
                        27, "PIPG_orig_precond_Dexplicit");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 27,
                        "PIPG_orig_precond_Dexplicit");
  }
  /* Call the function. */
  PIPG_orig_precond_Dexplicit_api(SD, prhs, nlhs, outputs);
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
  c_PIPG_orig_precond_DexplicitSt *d_PIPG_orig_precond_DexplicitSt = NULL;
  d_PIPG_orig_precond_DexplicitSt =
      (c_PIPG_orig_precond_DexplicitSt *)emlrtMxCalloc(
          (size_t)1, (size_t)1U * sizeof(c_PIPG_orig_precond_DexplicitSt));
  mexAtExit(&PIPG_orig_precond_Dexplicit_atexit);
  PIPG_orig_precond_Dexplicit_initialize();
  PIPG_orig_precond_Dexplicit_mexFunction(d_PIPG_orig_precond_DexplicitSt, nlhs,
                                          plhs, nrhs, prhs);
  PIPG_orig_precond_Dexplicit_terminate();
  emlrtMxFree(d_PIPG_orig_precond_DexplicitSt);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_PIPG_orig_precond_Dexplicit_mex.c) */
