/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_orig_precond_Dexplicit_types.h
 *
 * Code generation for function 'PIPG_orig_precond_Dexplicit'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"

/* Type Definitions */
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

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T
struct emxArray_int32_T {
  int32_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};
#endif /* struct_emxArray_int32_T */
#ifndef typedef_emxArray_int32_T
#define typedef_emxArray_int32_T
typedef struct emxArray_int32_T emxArray_int32_T;
#endif /* typedef_emxArray_int32_T */

#ifndef typedef_sparse
#define typedef_sparse
typedef struct {
  emxArray_real_T *d;
  emxArray_int32_T *colidx;
  emxArray_int32_T *rowidx;
} sparse;
#endif /* typedef_sparse */

#ifndef struct_emxArray_char_T_1x10
#define struct_emxArray_char_T_1x10
struct emxArray_char_T_1x10 {
  char_T data[10];
  int32_T size[2];
};
#endif /* struct_emxArray_char_T_1x10 */
#ifndef typedef_emxArray_char_T_1x10
#define typedef_emxArray_char_T_1x10
typedef struct emxArray_char_T_1x10 emxArray_char_T_1x10;
#endif /* typedef_emxArray_char_T_1x10 */

#ifndef typedef_rtString
#define typedef_rtString
typedef struct {
  emxArray_char_T_1x10 Value;
} rtString;
#endif /* typedef_rtString */

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
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_Ballz
#define typedef_Ballz
typedef struct {
  real_T indices[2];
  real_T r;
} Ballz;
#endif /* typedef_Ballz */

#ifndef typedef_Singletonz
#define typedef_Singletonz
typedef struct {
  real_T indices[3];
  real_T y[3];
} Singletonz;
#endif /* typedef_Singletonz */

#ifndef c_typedef_b_PIPG_orig_precond_D
#define c_typedef_b_PIPG_orig_precond_D
typedef struct {
  real_T S_z[1000000];
} b_PIPG_orig_precond_Dexplicit;
#endif /* c_typedef_b_PIPG_orig_precond_D */

#ifndef c_typedef_c_PIPG_orig_precond_D
#define c_typedef_c_PIPG_orig_precond_D
typedef struct {
  b_PIPG_orig_precond_Dexplicit f0;
} c_PIPG_orig_precond_DexplicitSt;
#endif /* c_typedef_c_PIPG_orig_precond_D */

/* End of code generation (PIPG_orig_precond_Dexplicit_types.h) */
