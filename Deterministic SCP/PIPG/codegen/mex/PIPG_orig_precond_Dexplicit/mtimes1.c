/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes1.c
 *
 * Code generation for function 'mtimes1'
 *
 */

/* Include files */
#include "mtimes1.h"
#include "PIPG_orig_precond_Dexplicit_data.h"
#include "PIPG_orig_precond_Dexplicit_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void sparse_mtimes(const emlrtStack *sp, const emxArray_real_T *a_d,
                   const emxArray_int32_T *a_colidx,
                   const emxArray_int32_T *a_rowidx, const real_T b[1000],
                   real_T c[3])
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  const real_T *a_d_data;
  const int32_T *a_colidx_data;
  const int32_T *a_rowidx_data;
  int32_T acol;
  int32_T ap;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  a_rowidx_data = a_rowidx->data;
  a_colidx_data = a_colidx->data;
  a_d_data = a_d->data;
  st.site = &q_emlrtRSI;
  c[0] = 0.0;
  c[1] = 0.0;
  c[2] = 0.0;
  if (a_colidx_data[a_colidx->size[0] - 1] - 1 != 0) {
    b_st.site = &r_emlrtRSI;
    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = 0.0;
    for (acol = 0; acol < 1000; acol++) {
      real_T bc;
      int32_T i;
      int32_T i1;
      int32_T nap;
      bc = b[acol];
      i = a_colidx_data[acol];
      i1 = a_colidx_data[acol + 1];
      nap = i1 - a_colidx_data[acol];
      if (nap >= 4) {
        int32_T apend1;
        apend1 = ((i1 - nap) + ((nap / 4) << 2)) - 1;
        c_st.site = &s_emlrtRSI;
        if ((a_colidx_data[acol] <= apend1) && (apend1 > 2147483643)) {
          d_st.site = &t_emlrtRSI;
          check_forloop_overflow_error(&d_st);
        }
        for (ap = i; ap <= apend1; ap += 4) {
          nap = a_rowidx_data[ap - 1] - 1;
          c[nap] += a_d_data[ap - 1] * bc;
          c[a_rowidx_data[ap] - 1] += a_d_data[ap] * bc;
          nap = a_rowidx_data[ap + 1] - 1;
          c[nap] += a_d_data[ap + 1] * bc;
          nap = a_rowidx_data[ap + 2] - 1;
          c[nap] += a_d_data[ap + 2] * bc;
        }
        i = apend1 + 1;
        for (ap = i; ap < i1; ap++) {
          nap = a_rowidx_data[ap - 1] - 1;
          c[nap] += a_d_data[ap - 1] * bc;
        }
      } else {
        for (ap = i; ap < i1; ap++) {
          nap = a_rowidx_data[ap - 1] - 1;
          c[nap] += a_d_data[ap - 1] * bc;
        }
      }
    }
  }
}

/* End of code generation (mtimes1.c) */
