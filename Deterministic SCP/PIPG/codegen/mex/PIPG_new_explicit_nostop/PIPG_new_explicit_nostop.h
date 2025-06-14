/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_new_explicit_nostop.h
 *
 * Code generation for function 'PIPG_new_explicit_nostop'
 *
 */

#pragma once

/* Include files */
#include "PIPG_new_explicit_nostop_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void PIPG_new_explicit_nostop(
    const emlrtStack *sp, const real_T qhat[1000], const real_T Hhat[3000],
    const real_T h[3], const real_T L[1000000], const real_T L_inv[1000000],
    real_T lambda, real_T sigma, real_T omega, real_T rho, real_T tol_abs,
    real_T tol_rel, real_T tol_infeas, real_T j_check, real_T j_max,
    const real_T z_ref[1000], const real_T what_ref[3], real_T z_star[1000],
    real_T what_star[3], struct0_T *sol_info);

/* End of code generation (PIPG_new_explicit_nostop.h) */
