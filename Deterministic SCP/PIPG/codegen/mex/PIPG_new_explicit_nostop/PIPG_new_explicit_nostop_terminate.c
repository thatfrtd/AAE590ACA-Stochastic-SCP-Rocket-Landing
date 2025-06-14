/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_new_explicit_nostop_terminate.c
 *
 * Code generation for function 'PIPG_new_explicit_nostop_terminate'
 *
 */

/* Include files */
#include "PIPG_new_explicit_nostop_terminate.h"
#include "PIPG_new_explicit_nostop_data.h"
#include "_coder_PIPG_new_explicit_nostop_mex.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void emlrtExitTimeCleanupDtorFcn(const void *r);

/* Function Definitions */
static void emlrtExitTimeCleanupDtorFcn(const void *r)
{
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void PIPG_new_explicit_nostop_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtPushHeapReferenceStackR2021a(
      &st, false, NULL, (void *)&emlrtExitTimeCleanupDtorFcn, NULL, NULL, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtProfilerUnregisterMEXFcn((char_T *)c_PIPG_new_explicit_nostop_comp,
                                isMexOutdated);
  emlrtProfilerUnregisterMEXFcn((char_T *)ball_project_complete_name,
                                isMexOutdated);
  emlrtProfilerUnregisterMEXFcn((char_T *)singleton_project_complete_name,
                                isMexOutdated);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void PIPG_new_explicit_nostop_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (PIPG_new_explicit_nostop_terminate.c) */
