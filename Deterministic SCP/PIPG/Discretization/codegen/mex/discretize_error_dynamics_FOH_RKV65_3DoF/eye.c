/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eye.c
 *
 * Code generation for function 'eye'
 *
 */

/* Include files */
#include "eye.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void eye(real_T b_I[49])
{
  int32_T k;
  memset(&b_I[0], 0, 49U * sizeof(real_T));
  for (k = 0; k < 7; k++) {
    b_I[k + 7 * k] = 1.0;
  }
}

/* End of code generation (eye.c) */
