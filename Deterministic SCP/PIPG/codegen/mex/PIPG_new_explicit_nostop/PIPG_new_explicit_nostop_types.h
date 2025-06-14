/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_new_explicit_nostop_types.h
 *
 * Code generation for function 'PIPG_new_explicit_nostop'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef typedef_rtString
#define typedef_rtString
typedef struct {
  char_T Value[8];
} rtString;
#endif /* typedef_rtString */

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T {
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};
#endif /* struct_emxArray_real_T */
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /* typedef_emxArray_real_T */

#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  boolean_T primal_infeasible;
  boolean_T dual_infeasible;
  rtString solution_status;
  real_T iterations;
  emxArray_real_T *zhat_inf_dj;
  emxArray_real_T *what_inf_dj;
  emxArray_real_T *zhis;
  real_T time;
} struct0_T;
#endif /* typedef_struct0_T */

/* End of code generation (PIPG_new_explicit_nostop_types.h) */
