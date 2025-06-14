/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * discretize_error_dynamics_FOH_RKV65_3DoF.h
 *
 * Code generation for function 'discretize_error_dynamics_FOH_RKV65_3DoF'
 *
 */

#pragma once

/* Include files */
#include "discretize_error_dynamics_FOH_RKV65_3DoF_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void discretize_error_dynamics_FOH_RKV65_3DoF(
    const emlrtStack *sp, real_T N, const real_T tspan[2],
    const real_T x_ref[105], const real_T u_ref[30], real_T s_ref, real_T N_sub,
    real_T L, real_T b_I, real_T alpha, emxArray_real_T *A_k,
    emxArray_real_T *B_k_plus, emxArray_real_T *B_k_minus, emxArray_real_T *S_k,
    emxArray_real_T *d_k, emxArray_real_T *Delta);

/* End of code generation (discretize_error_dynamics_FOH_RKV65_3DoF.h) */
