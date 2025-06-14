/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * custom_power_iteration_noops.h
 *
 * Code generation for function 'custom_power_iteration_noops'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
real_T custom_power_iteration_noops(
    const emlrtStack *sp, const real_T Ahat_1[686], const real_T Ahat_2[686],
    const real_T Bhat_minus[196], const real_T Bhat_plus[196],
    const real_T Shat[98], real_T lx1_inv, real_T lx2_inv, real_T x[105],
    real_T xi[105], real_T u[30], real_T s);

/* End of code generation (custom_power_iteration_noops.h) */
