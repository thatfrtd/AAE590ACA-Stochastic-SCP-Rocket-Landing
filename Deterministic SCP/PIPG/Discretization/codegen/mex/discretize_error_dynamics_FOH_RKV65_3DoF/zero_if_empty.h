/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * zero_if_empty.h
 *
 * Code generation for function 'zero_if_empty'
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
void zero_if_empty(const emlrtStack *sp, const emxArray_real_T *v,
                   emxArray_real_T *b_v);

/* End of code generation (zero_if_empty.h) */
