/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse.c
 *
 * Code generation for function 'sparse'
 *
 */

/* Include files */
#include "sparse.h"
#include "PIPG_orig_precond_Dexplicit_data.h"
#include "PIPG_orig_precond_Dexplicit_emxutil.h"
#include "PIPG_orig_precond_Dexplicit_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"
#include <emmintrin.h>

/* Variable Definitions */
static emlrtRSInfo u_emlrtRSI = {
    395,                 /* lineNo */
    "sparse/ctranspose", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo v_emlrtRSI = {
    17,                    /* lineNo */
    "sparse/locTranspose", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\locTranspose.m" /* pathName */
};

static emlrtRSInfo w_emlrtRSI = {
    176,             /* lineNo */
    "sparse/sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo ab_emlrtRSI = {
    302,            /* lineNo */
    "sparse/times", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo bb_emlrtRSI = {
    194,            /* lineNo */
    "sparse/binOp", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\binOp.m" /* pathName */
};

static emlrtRSInfo cb_emlrtRSI = {
    204,            /* lineNo */
    "sparse/binOp", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\binOp.m" /* pathName */
};

static emlrtRSInfo db_emlrtRSI = {
    206,            /* lineNo */
    "sparse/binOp", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\binOp.m" /* pathName */
};

static emlrtRSInfo eb_emlrtRSI = {
    18,      /* lineNo */
    "spfun", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\sparfun\\spfun.m" /* pathName
                                                                          */
};

static emlrtRSInfo fb_emlrtRSI = {
    462,                /* lineNo */
    "sparse/spfunImpl", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo gb_emlrtRSI = {
    465,                /* lineNo */
    "sparse/spfunImpl", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo hb_emlrtRSI = {
    468,                /* lineNo */
    "sparse/spfunImpl", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo ib_emlrtRSI = {
    1476,                 /* lineNo */
    "sparse/spallocLike", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo kb_emlrtRSI = {
    14,              /* lineNo */
    "sparse/fillIn", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\fillIn.m" /* pathName */
};

static emlrtRSInfo lb_emlrtRSI = {
    437,           /* lineNo */
    "scalarBinOp", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\binOp.m" /* pathName */
};

static emlrtRSInfo mb_emlrtRSI = {
    13,       /* lineNo */
    "sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\sparfun\\sparse.m" /* pathName
                                                                           */
};

static emlrtRTEInfo d_emlrtRTEI = {
    178,             /* lineNo */
    39,              /* colNo */
    "sparse/sparse", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo e_emlrtRTEI = {
    1623,              /* lineNo */
    9,                 /* colNo */
    "assertValidSize", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo f_emlrtRTEI = {
    1626,              /* lineNo */
    31,                /* colNo */
    "assertValidSize", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo g_emlrtRTEI = {
    460,                /* lineNo */
    34,                 /* colNo */
    "sparse/spfunImpl", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo m_emlrtRTEI = {
    395,      /* lineNo */
    13,       /* colNo */
    "sparse", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo n_emlrtRTEI = {
    203,     /* lineNo */
    9,       /* colNo */
    "binOp", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\binOp.m" /* pName */
};

static emlrtRTEInfo o_emlrtRTEI = {
    459,      /* lineNo */
    12,       /* colNo */
    "sparse", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo p_emlrtRTEI = {
    302,      /* lineNo */
    13,       /* colNo */
    "sparse", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo q_emlrtRTEI = {
    125,     /* lineNo */
    5,       /* colNo */
    "binOp", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\binOp.m" /* pName */
};

/* Function Definitions */
void sparse_ctranspose(const emlrtStack *sp, const emxArray_real_T *this_d,
                       const emxArray_int32_T *this_colidx,
                       const emxArray_int32_T *this_rowidx,
                       emxArray_real_T *y_d, emxArray_int32_T *y_colidx,
                       emxArray_int32_T *y_rowidx)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  const real_T *this_d_data;
  real_T *y_d_data;
  int32_T counts[3];
  const int32_T *this_colidx_data;
  const int32_T *this_rowidx_data;
  int32_T c;
  int32_T numalloc;
  int32_T outridx;
  int32_T y_tmp_tmp;
  int32_T *y_colidx_data;
  int32_T *y_rowidx_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  this_rowidx_data = this_rowidx->data;
  this_colidx_data = this_colidx->data;
  this_d_data = this_d->data;
  st.site = &u_emlrtRSI;
  y_tmp_tmp = this_colidx_data[this_colidx->size[0] - 1];
  b_st.site = &v_emlrtRSI;
  c_st.site = &w_emlrtRSI;
  if (y_tmp_tmp - 1 < 0) {
    emlrtErrorWithMessageIdR2018a(&c_st, &e_emlrtRTEI,
                                  "Coder:toolbox:SparseNegativeSize",
                                  "Coder:toolbox:SparseNegativeSize", 0);
  }
  if (y_tmp_tmp - 1 < 0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
                                  "Coder:toolbox:SparseNzmaxTooSmall",
                                  "Coder:toolbox:SparseNzmaxTooSmall", 0);
  }
  if (y_tmp_tmp - 1 >= 1) {
    numalloc = y_tmp_tmp - 2;
  } else {
    numalloc = 0;
  }
  outridx = y_d->size[0];
  y_d->size[0] = numalloc + 1;
  emxEnsureCapacity_real_T(&b_st, y_d, outridx, &m_emlrtRTEI);
  y_d_data = y_d->data;
  for (outridx = 0; outridx <= numalloc; outridx++) {
    y_d_data[outridx] = 0.0;
  }
  outridx = y_rowidx->size[0];
  y_rowidx->size[0] = numalloc + 1;
  emxEnsureCapacity_int32_T(&b_st, y_rowidx, outridx, &m_emlrtRTEI);
  y_rowidx_data = y_rowidx->data;
  for (outridx = 0; outridx <= numalloc; outridx++) {
    y_rowidx_data[outridx] = 0;
  }
  outridx = y_colidx->size[0];
  y_colidx->size[0] = 4;
  emxEnsureCapacity_int32_T(&st, y_colidx, outridx, &m_emlrtRTEI);
  y_colidx_data = y_colidx->data;
  for (outridx = 0; outridx < 4; outridx++) {
    y_colidx_data[outridx] = 0;
  }
  for (numalloc = 0; numalloc <= y_tmp_tmp - 2; numalloc++) {
    y_colidx_data[this_rowidx_data[numalloc]]++;
  }
  y_colidx_data[0] = 1;
  y_colidx_data[1] += y_colidx_data[0];
  counts[0] = 0;
  y_colidx_data[2] += y_colidx_data[1];
  counts[1] = 0;
  y_colidx_data[3] += y_colidx_data[2];
  counts[2] = 0;
  for (c = 0; c < 1000; c++) {
    for (numalloc = this_colidx_data[c] - 1;
         numalloc + 1 < this_colidx_data[c + 1]; numalloc++) {
      y_tmp_tmp = counts[this_rowidx_data[numalloc] - 1];
      outridx = (y_tmp_tmp + y_colidx_data[this_rowidx_data[numalloc] - 1]) - 1;
      y_d_data[outridx] = this_d_data[numalloc];
      y_rowidx_data[outridx] = c + 1;
      counts[this_rowidx_data[numalloc] - 1] = y_tmp_tmp + 1;
    }
  }
}

void sparse_times(const emlrtStack *sp, real_T a, const emxArray_real_T *b_d,
                  const emxArray_int32_T *b_colidx,
                  const emxArray_int32_T *b_rowidx, emxArray_real_T *s_d,
                  emxArray_int32_T *s_colidx, emxArray_int32_T *s_rowidx)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  emxArray_real_T *S;
  const real_T *b_d_data;
  real_T *S_data;
  real_T *s_d_data;
  const int32_T *b_colidx_data;
  const int32_T *b_rowidx_data;
  int32_T currRowIdx;
  int32_T i;
  int32_T idx;
  int32_T nzs_tmp_tmp;
  int32_T *s_colidx_data;
  int32_T *s_rowidx_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  b_rowidx_data = b_rowidx->data;
  b_colidx_data = b_colidx->data;
  b_d_data = b_d->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  st.site = &ab_emlrtRSI;
  if (a * 0.0 == 0.0) {
    int32_T loop_ub;
    int32_T numalloc;
    b_st.site = &bb_emlrtRSI;
    c_st.site = &eb_emlrtRSI;
    nzs_tmp_tmp = b_colidx_data[b_colidx->size[0] - 1];
    if (nzs_tmp_tmp - 1 < 1) {
      loop_ub = 0;
    } else {
      loop_ub = nzs_tmp_tmp - 1;
    }
    emxInit_real_T(&c_st, &S, 1, &o_emlrtRTEI);
    i = S->size[0];
    S->size[0] = loop_ub;
    emxEnsureCapacity_real_T(&c_st, S, i, &o_emlrtRTEI);
    S_data = S->data;
    numalloc = (loop_ub / 2) << 1;
    currRowIdx = numalloc - 2;
    for (i = 0; i <= currRowIdx; i += 2) {
      _mm_storeu_pd(&S_data[i],
                    _mm_mul_pd(_mm_set1_pd(a), _mm_loadu_pd(&b_d_data[i])));
    }
    for (i = numalloc; i < loop_ub; i++) {
      S_data[i] = a * b_d_data[i];
    }
    if (S->size[0] != nzs_tmp_tmp - 1) {
      emlrtErrorWithMessageIdR2018a(&c_st, &g_emlrtRTEI, "MATLAB:samelen",
                                    "MATLAB:samelen", 0);
    }
    d_st.site = &fb_emlrtRSI;
    e_st.site = &ib_emlrtRSI;
    f_st.site = &w_emlrtRSI;
    if (nzs_tmp_tmp - 1 < 0) {
      emlrtErrorWithMessageIdR2018a(&f_st, &e_emlrtRTEI,
                                    "Coder:toolbox:SparseNegativeSize",
                                    "Coder:toolbox:SparseNegativeSize", 0);
    }
    if (nzs_tmp_tmp - 1 >= MAX_int32_T) {
      emlrtErrorWithMessageIdR2018a(
          &f_st, &f_emlrtRTEI, "Coder:toolbox:SparseMaxSize",
          "Coder:toolbox:SparseMaxSize", 2, 12, MAX_int32_T);
    }
    if (nzs_tmp_tmp - 1 < 0) {
      emlrtErrorWithMessageIdR2018a(&e_st, &d_emlrtRTEI,
                                    "Coder:toolbox:SparseNzmaxTooSmall",
                                    "Coder:toolbox:SparseNzmaxTooSmall", 0);
    }
    if (nzs_tmp_tmp - 1 >= 1) {
      numalloc = nzs_tmp_tmp - 2;
    } else {
      numalloc = 0;
    }
    i = s_d->size[0];
    s_d->size[0] = numalloc + 1;
    emxEnsureCapacity_real_T(&e_st, s_d, i, &p_emlrtRTEI);
    s_d_data = s_d->data;
    i = s_rowidx->size[0];
    s_rowidx->size[0] = numalloc + 1;
    emxEnsureCapacity_int32_T(&e_st, s_rowidx, i, &p_emlrtRTEI);
    s_rowidx_data = s_rowidx->data;
    for (i = 0; i <= numalloc; i++) {
      s_d_data[i] = 0.0;
      s_rowidx_data[i] = 0;
    }
    if (nzs_tmp_tmp - 1 < 1) {
      loop_ub = 1;
    } else {
      loop_ub = nzs_tmp_tmp;
    }
    for (i = 0; i <= loop_ub - 2; i++) {
      s_rowidx_data[i] = b_rowidx_data[i];
    }
    loop_ub = b_colidx->size[0];
    i = s_colidx->size[0];
    s_colidx->size[0] = b_colidx->size[0];
    emxEnsureCapacity_int32_T(&c_st, s_colidx, i, &p_emlrtRTEI);
    s_colidx_data = s_colidx->data;
    for (i = 0; i < loop_ub; i++) {
      s_colidx_data[i] = b_colidx_data[i];
    }
    d_st.site = &gb_emlrtRSI;
    if (nzs_tmp_tmp - 1 > 2147483646) {
      e_st.site = &t_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }
    for (currRowIdx = 0; currRowIdx <= nzs_tmp_tmp - 2; currRowIdx++) {
      s_d_data[currRowIdx] = S_data[currRowIdx];
    }
    emxFree_real_T(&c_st, &S);
    d_st.site = &hb_emlrtRSI;
    idx = 1;
    e_st.site = &kb_emlrtRSI;
    for (nzs_tmp_tmp = 0; nzs_tmp_tmp <= loop_ub - 2; nzs_tmp_tmp++) {
      numalloc = s_colidx_data[nzs_tmp_tmp];
      s_colidx_data[nzs_tmp_tmp] = idx;
      while (numalloc < s_colidx_data[nzs_tmp_tmp + 1]) {
        real_T val;
        currRowIdx = s_rowidx_data[numalloc - 1];
        val = s_d_data[numalloc - 1];
        numalloc++;
        if (val != 0.0) {
          s_d_data[idx - 1] = val;
          s_rowidx_data[idx - 1] = currRowIdx;
          idx++;
        }
      }
    }
    s_colidx_data[b_colidx->size[0] - 1] = idx;
  } else {
    int32_T numalloc;
    emxInit_real_T(&st, &S, 2, &q_emlrtRTEI);
    i = S->size[0] * S->size[1];
    S->size[0] = 3;
    S->size[1] = 1000;
    emxEnsureCapacity_real_T(&st, S, i, &n_emlrtRTEI);
    S_data = S->data;
    for (i = 0; i < 3000; i++) {
      S_data[i] = rtNaN;
    }
    b_st.site = &cb_emlrtRSI;
    c_st.site = &lb_emlrtRSI;
    for (currRowIdx = 0; currRowIdx < 1000; currRowIdx++) {
      i = b_colidx_data[currRowIdx];
      numalloc = b_colidx_data[currRowIdx + 1];
      for (idx = i; idx < numalloc; idx++) {
        S_data[(b_rowidx_data[idx - 1] + S->size[0] * currRowIdx) - 1] =
            a * b_d_data[idx - 1];
      }
    }
    b_st.site = &db_emlrtRSI;
    c_st.site = &mb_emlrtRSI;
    numalloc = 0;
    for (currRowIdx = 0; currRowIdx < 3000; currRowIdx++) {
      if (S_data[currRowIdx % 3 + S->size[0] * (currRowIdx / 3)] != 0.0) {
        numalloc++;
      }
    }
    numalloc = muIntScalarMax_sint32(numalloc, 1);
    i = s_d->size[0];
    s_d->size[0] = numalloc;
    emxEnsureCapacity_real_T(&c_st, s_d, i, &p_emlrtRTEI);
    s_d_data = s_d->data;
    for (i = 0; i < numalloc; i++) {
      s_d_data[i] = 0.0;
    }
    i = s_colidx->size[0];
    s_colidx->size[0] = 1001;
    emxEnsureCapacity_int32_T(&c_st, s_colidx, i, &p_emlrtRTEI);
    s_colidx_data = s_colidx->data;
    for (i = 0; i < 1001; i++) {
      s_colidx_data[i] = 0;
    }
    s_colidx_data[0] = 1;
    i = s_rowidx->size[0];
    s_rowidx->size[0] = numalloc;
    emxEnsureCapacity_int32_T(&c_st, s_rowidx, i, &p_emlrtRTEI);
    s_rowidx_data = s_rowidx->data;
    for (i = 0; i < numalloc; i++) {
      s_rowidx_data[i] = 0;
    }
    s_rowidx_data[0] = 1;
    numalloc = 0;
    for (currRowIdx = 0; currRowIdx < 1000; currRowIdx++) {
      real_T val;
      val = S_data[S->size[0] * currRowIdx];
      if (val != 0.0) {
        s_rowidx_data[numalloc] = 1;
        s_d_data[numalloc] = val;
        numalloc++;
      }
      val = S_data[S->size[0] * currRowIdx + 1];
      if (val != 0.0) {
        s_rowidx_data[numalloc] = 2;
        s_d_data[numalloc] = val;
        numalloc++;
      }
      val = S_data[S->size[0] * currRowIdx + 2];
      if (val != 0.0) {
        s_rowidx_data[numalloc] = 3;
        s_d_data[numalloc] = val;
        numalloc++;
      }
      s_colidx_data[currRowIdx + 1] = numalloc + 1;
    }
    emxFree_real_T(&c_st, &S);
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (sparse.c) */
