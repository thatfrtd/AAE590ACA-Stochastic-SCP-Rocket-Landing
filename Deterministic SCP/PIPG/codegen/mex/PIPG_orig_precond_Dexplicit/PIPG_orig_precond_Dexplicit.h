/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_orig_precond_Dexplicit.h
 *
 * Code generation for function 'PIPG_orig_precond_Dexplicit'
 *
 */

#pragma once

/* Include files */
#include "PIPG_orig_precond_Dexplicit_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void PIPG_orig_precond_Dexplicit(
    c_PIPG_orig_precond_DexplicitSt *SD, const emlrtStack *sp,
    const real_T qhat[1000], const sparse Hhat, const real_T h[3],
    const Ballz *Dball, const Singletonz *Dsingleton, const real_T L[1000000],
    const real_T L_inv[1000000], real_T lambda, real_T sigma, real_T rho,
    real_T tol_abs, real_T tol_rel, real_T tol_infeas, real_T j_check,
    real_T j_restart, real_T j_max, const real_T z_ref[1000],
    const real_T what_ref[3], real_T S_z, real_T c_z, real_T z_star[1000],
    real_T what_star[3], struct0_T *sol_info);

/* End of code generation (PIPG_orig_precond_Dexplicit.h) */
