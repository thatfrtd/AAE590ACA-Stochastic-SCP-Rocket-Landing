/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * stopping.h
 *
 * Code generation for function 'stopping'
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
boolean_T stopping(const real_T zhat_jp1[1000], const real_T what_jp1[3],
                   const real_T zhat_j[1000], const real_T what_j[3],
                   real_T tol_abs, real_T tol_rel, real_T *zhat_inf_dj,
                   real_T *what_inf_dj);

/* End of code generation (stopping.h) */
