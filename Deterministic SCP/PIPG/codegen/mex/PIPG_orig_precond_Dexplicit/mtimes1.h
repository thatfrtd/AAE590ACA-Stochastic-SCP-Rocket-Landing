/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes1.h
 *
 * Code generation for function 'mtimes1'
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
void sparse_mtimes(const emlrtStack *sp, const emxArray_real_T *a_d,
                   const emxArray_int32_T *a_colidx,
                   const emxArray_int32_T *a_rowidx, const real_T b[1000],
                   real_T c[3]);

/* End of code generation (mtimes1.h) */
