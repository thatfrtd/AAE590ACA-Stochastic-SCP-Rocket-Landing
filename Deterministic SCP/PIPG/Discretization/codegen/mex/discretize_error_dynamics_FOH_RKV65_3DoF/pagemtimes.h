/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * pagemtimes.h
 *
 * Code generation for function 'pagemtimes'
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
void b_pagemtimes(const emlrtStack *sp, const emxArray_real_T *varargin_1,
                  const emxArray_real_T *varargin_2, emxArray_real_T *Z);

void c_pagemtimes(const emlrtStack *sp, const emxArray_real_T *varargin_1,
                  const emxArray_real_T *varargin_2, emxArray_real_T *Z);

void d_pagemtimes(const emlrtStack *sp, const emxArray_real_T *varargin_1,
                  const real_T varargin_2_data[],
                  const int32_T varargin_2_size[3], emxArray_real_T *Z);

void pagemtimes(const emlrtStack *sp, const emxArray_real_T *varargin_1,
                const emxArray_real_T *varargin_2, emxArray_real_T *Z);

/* End of code generation (pagemtimes.h) */
