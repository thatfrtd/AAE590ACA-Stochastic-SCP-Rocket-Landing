/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_orig_precond_Dexplicit_initialize.c
 *
 * Code generation for function 'PIPG_orig_precond_Dexplicit_initialize'
 *
 */

/* Include files */
#include "PIPG_orig_precond_Dexplicit_initialize.h"
#include "PIPG_orig_precond_Dexplicit_data.h"
#include "_coder_PIPG_orig_precond_Dexplicit_mex.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void PIPG_orig_precond_Dexplicit_once(void);

/* Function Definitions */
static void PIPG_orig_precond_Dexplicit_once(void)
{
  mex_InitInfAndNan();
}

void PIPG_orig_precond_Dexplicit_initialize(void)
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
    PIPG_orig_precond_Dexplicit_once();
  }
}

/* End of code generation (PIPG_orig_precond_Dexplicit_initialize.c) */
