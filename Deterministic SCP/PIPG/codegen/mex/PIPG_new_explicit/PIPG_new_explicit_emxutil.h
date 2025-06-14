/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_new_explicit_emxutil.h
 *
 * Code generation for function 'PIPG_new_explicit_emxutil'
 *
 */

#pragma once

/* Include files */
#include "PIPG_new_explicit_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void emxEnsureCapacity_real_T(const emlrtStack *sp, emxArray_real_T *emxArray,
                              int32_T oldNumel,
                              const emlrtRTEInfo *srcLocation);

void emxFreeStruct_struct0_T(const emlrtStack *sp, struct0_T *pStruct);

void emxFree_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray);

void emxInitStruct_rtString(rtString *pStruct);

void emxInitStruct_struct0_T(const emlrtStack *sp, struct0_T *pStruct,
                             const emlrtRTEInfo *srcLocation);

void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
                    int32_T numDimensions, const emlrtRTEInfo *srcLocation);

/* End of code generation (PIPG_new_explicit_emxutil.h) */
