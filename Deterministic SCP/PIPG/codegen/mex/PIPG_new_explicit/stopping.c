/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * stopping.c
 *
 * Code generation for function 'stopping'
 *
 */

/* Include files */
#include "stopping.h"
#include "PIPG_new_explicit_data.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
boolean_T stopping(const real_T zhat_jp1[1000], const real_T what_jp1[3],
                   const real_T zhat_j[1000], const real_T what_j[3],
                   real_T tol_abs, real_T tol_rel, real_T *zhat_inf_dj,
                   real_T *what_inf_dj)
{
  real_T absx;
  real_T what_inf_j;
  real_T what_inf_jp1;
  real_T zhat_inf_j;
  real_T zhat_inf_jp1;
  int32_T k;
  boolean_T terminate;
  emlrtMEXProfilingFunctionEntry((char_T *)stopping_complete_name,
                                 isMexOutdated);
  /* STOPPING Summary of this function goes here */
  /*    Detailed explanation goes here */
  emlrtMEXProfilingStatement(1, isMexOutdated);
  zhat_inf_jp1 = 0.0;
  for (k = 0; k < 1000; k++) {
    absx = muDoubleScalarAbs(zhat_jp1[k]);
    if (muDoubleScalarIsNaN(absx) || (absx > zhat_inf_jp1)) {
      zhat_inf_jp1 = absx;
    }
  }
  emlrtMEXProfilingStatement(2, isMexOutdated);
  zhat_inf_j = 0.0;
  for (k = 0; k < 1000; k++) {
    absx = muDoubleScalarAbs(zhat_j[k]);
    if (muDoubleScalarIsNaN(absx) || (absx > zhat_inf_j)) {
      zhat_inf_j = absx;
    }
  }
  emlrtMEXProfilingStatement(3, isMexOutdated);
  *zhat_inf_dj = 0.0;
  for (k = 0; k < 1000; k++) {
    absx = muDoubleScalarAbs(zhat_jp1[k] - zhat_j[k]);
    if (muDoubleScalarIsNaN(absx) || (absx > *zhat_inf_dj)) {
      *zhat_inf_dj = absx;
    }
  }
  emlrtMEXProfilingStatement(4, isMexOutdated);
  what_inf_jp1 = 0.0;
  for (k = 0; k < 3; k++) {
    absx = muDoubleScalarAbs(what_jp1[k]);
    if (muDoubleScalarIsNaN(absx) || (absx > what_inf_jp1)) {
      what_inf_jp1 = absx;
    }
  }
  emlrtMEXProfilingStatement(5, isMexOutdated);
  what_inf_j = 0.0;
  for (k = 0; k < 3; k++) {
    absx = muDoubleScalarAbs(what_j[k]);
    if (muDoubleScalarIsNaN(absx) || (absx > what_inf_j)) {
      what_inf_j = absx;
    }
  }
  emlrtMEXProfilingStatement(6, isMexOutdated);
  *what_inf_dj = 0.0;
  for (k = 0; k < 3; k++) {
    absx = muDoubleScalarAbs(what_jp1[k] - what_j[k]);
    if (muDoubleScalarIsNaN(absx) || (absx > *what_inf_dj)) {
      *what_inf_dj = absx;
    }
  }
  emlrtMEXProfilingStatement(7, isMexOutdated);
  if ((*zhat_inf_dj <=
       tol_abs + tol_rel * muDoubleScalarMax(zhat_inf_jp1, zhat_inf_j)) &&
      (*what_inf_dj <=
       tol_abs + tol_rel * muDoubleScalarMax(what_inf_jp1, what_inf_j))) {
    emlrtMEXProfilingStatement(8, isMexOutdated);
    terminate = true;
    emlrtMEXProfilingStatement(11, isMexOutdated);
  } else {
    emlrtMEXProfilingStatement(9, isMexOutdated);
    emlrtMEXProfilingStatement(10, isMexOutdated);
    terminate = false;
    emlrtMEXProfilingStatement(11, isMexOutdated);
  }
  emlrtMEXProfilingStatement(12, isMexOutdated);
  emlrtMEXProfilingFunctionExit(isMexOutdated);
  return terminate;
}

/* End of code generation (stopping.c) */
