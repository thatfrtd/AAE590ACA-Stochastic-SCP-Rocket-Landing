/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * linspace.c
 *
 * Code generation for function 'linspace'
 *
 */

/* Include files */
#include "linspace.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_emxutil.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRTEInfo f_emlrtRTEI = {
    33,         /* lineNo */
    37,         /* colNo */
    "linspace", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\elmat\\linspace.m" /* pName
                                                                           */
};

static emlrtRTEInfo uc_emlrtRTEI = {
    49,         /* lineNo */
    20,         /* colNo */
    "linspace", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\elmat\\linspace.m" /* pName
                                                                           */
};

/* Function Definitions */
void linspace(const emlrtStack *sp, real_T d1, real_T d2, real_T n,
              emxArray_real_T *y)
{
  real_T *y_data;
  int32_T k;
  if (!(n >= 0.0)) {
    if (muDoubleScalarIsNaN(n)) {
      emlrtErrorWithMessageIdR2018a(sp, &f_emlrtRTEI,
                                    "Coder:toolbox:MustNotBeNaN",
                                    "Coder:toolbox:MustNotBeNaN", 3, 4, 1, "N");
    }
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    real_T delta1;
    int32_T i;
    delta1 = muDoubleScalarFloor(n);
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    i = (int32_T)muDoubleScalarFloor(n);
    y->size[1] = (int32_T)delta1;
    emxEnsureCapacity_real_T(sp, y, k, &uc_emlrtRTEI);
    y_data = y->data;
    if ((int32_T)delta1 >= 1) {
      y_data[(int32_T)delta1 - 1] = d2;
      if (y->size[1] >= 2) {
        y_data[0] = d1;
        if (y->size[1] >= 3) {
          if (d1 == -d2) {
            delta1 = d2 / ((real_T)y->size[1] - 1.0);
            for (k = 2; k < i; k++) {
              y_data[k - 1] = (real_T)(((k << 1) - y->size[1]) - 1) * delta1;
            }
            if (((uint32_T)y->size[1] & 1U) == 1U) {
              y_data[y->size[1] >> 1] = 0.0;
            }
          } else if (((d1 < 0.0) != (d2 < 0.0)) &&
                     ((muDoubleScalarAbs(d1) > 8.9884656743115785E+307) ||
                      (muDoubleScalarAbs(d2) > 8.9884656743115785E+307))) {
            real_T delta2;
            delta1 = d1 / ((real_T)y->size[1] - 1.0);
            delta2 = d2 / ((real_T)y->size[1] - 1.0);
            for (k = 0; k <= i - 3; k++) {
              y_data[k + 1] = (d1 + delta2 * ((real_T)k + 1.0)) -
                              delta1 * ((real_T)k + 1.0);
            }
          } else {
            delta1 = (d2 - d1) / ((real_T)y->size[1] - 1.0);
            for (k = 0; k <= i - 3; k++) {
              y_data[k + 1] = d1 + ((real_T)k + 1.0) * delta1;
            }
          }
        }
      }
    }
  }
}

/* End of code generation (linspace.c) */
