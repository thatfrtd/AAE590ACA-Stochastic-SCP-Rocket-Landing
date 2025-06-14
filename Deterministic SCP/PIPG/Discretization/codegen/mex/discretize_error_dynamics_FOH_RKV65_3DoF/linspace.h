/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * linspace.h
 *
 * Code generation for function 'linspace'
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
void linspace(const emlrtStack *sp, real_T d1, real_T d2, real_T n,
              emxArray_real_T *y);

/* End of code generation (linspace.h) */
