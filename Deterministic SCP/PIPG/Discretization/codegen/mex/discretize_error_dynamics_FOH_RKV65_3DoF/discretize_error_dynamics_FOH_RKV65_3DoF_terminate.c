/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * discretize_error_dynamics_FOH_RKV65_3DoF_terminate.c
 *
 * Code generation for function
 * 'discretize_error_dynamics_FOH_RKV65_3DoF_terminate'
 *
 */

/* Include files */
#include "discretize_error_dynamics_FOH_RKV65_3DoF_terminate.h"
#include "_coder_discretize_error_dynamics_FOH_RKV65_3DoF_mex.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void emlrtExitTimeCleanupDtorFcn(const void *r);

/* Function Definitions */
static void emlrtExitTimeCleanupDtorFcn(const void *r)
{
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void discretize_error_dynamics_FOH_RKV65_3DoF_atexit(void)
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
  emlrtProfilerUnregisterMEXFcn((char_T *)d_discretize_error_dynamics_FOH,
                                isMexOutdated);
  emlrtProfilerUnregisterMEXFcn((char_T *)STM_diff_eq_FOH_complete_name,
                                isMexOutdated);
  emlrtProfilerUnregisterMEXFcn((char_T *)A_3DoF_complete_name, isMexOutdated);
  emlrtProfilerUnregisterMEXFcn((char_T *)B_3DoF_complete_name, isMexOutdated);
  emlrtProfilerUnregisterMEXFcn((char_T *)c_SymDynamics3DoF_mass_polar_co,
                                isMexOutdated);
  emlrtProfilerUnregisterMEXFcn((char_T *)S_3DoF_complete_name, isMexOutdated);
  emlrtProfilerUnregisterMEXFcn((char_T *)zero_if_empty_complete_name,
                                isMexOutdated);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void discretize_error_dynamics_FOH_RKV65_3DoF_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (discretize_error_dynamics_FOH_RKV65_3DoF_terminate.c)
 */
