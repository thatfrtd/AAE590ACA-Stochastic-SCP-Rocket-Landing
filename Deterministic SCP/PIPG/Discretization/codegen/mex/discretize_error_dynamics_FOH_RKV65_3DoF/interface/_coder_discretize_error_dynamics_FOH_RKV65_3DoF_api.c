/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_discretize_error_dynamics_FOH_RKV65_3DoF_api.c
 *
 * Code generation for function
 * '_coder_discretize_error_dynamics_FOH_RKV65_3DoF_api'
 *
 */

/* Include files */
#include "_coder_discretize_error_dynamics_FOH_RKV65_3DoF_api.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_data.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_emxutil.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo id_emlrtRTEI = {
    1,                                                     /* lineNo */
    1,                                                     /* colNo */
    "_coder_discretize_error_dynamics_FOH_RKV65_3DoF_api", /* fName */
    ""                                                     /* pName */
};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static const mxArray *b_emlrt_marshallOut(emxArray_real_T *u);

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[2];

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[2];

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[105];

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier);

static const mxArray *emlrt_marshallOut(emxArray_real_T *u);

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[105];

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[30];

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[30];

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[2];

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[105];

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[30];

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *b_emlrt_marshallOut(emxArray_real_T *u)
{
  static const int32_T iv[2] = {0, 0};
  const mxArray *m;
  const mxArray *y;
  real_T *u_data;
  u_data = u->data;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u_data[0]);
  emlrtSetDimensions((mxArray *)m, &u->size[0], 2);
  u->canFreeData = false;
  emlrtAssign(&y, m);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[2]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[2];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[2]
{
  real_T(*y)[2];
  y = j_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[105]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[105];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *emlrt_marshallOut(emxArray_real_T *u)
{
  static const int32_T iv[3] = {0, 0, 0};
  const mxArray *m;
  const mxArray *y;
  real_T *u_data;
  u_data = u->data;
  y = NULL;
  m = emlrtCreateNumericArray(3, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u_data[0]);
  emlrtSetDimensions((mxArray *)m, &u->size[0], 3);
  u->canFreeData = false;
  emlrtAssign(&y, m);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[105]
{
  real_T(*y)[105];
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[30]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[30];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[30]
{
  real_T(*y)[30];
  y = l_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[2]
{
  static const int32_T dims[2] = {1, 2};
  real_T(*ret)[2];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[2])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[105]
{
  static const int32_T dims[2] = {7, 15};
  real_T(*ret)[105];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[105])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[30]
{
  static const int32_T dims[2] = {2, 15};
  real_T(*ret)[30];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[30])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void e_discretize_error_dynamics_FOH(const mxArray *const prhs[9], int32_T nlhs,
                                     const mxArray *plhs[6])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  emxArray_real_T *A_k;
  emxArray_real_T *B_k_minus;
  emxArray_real_T *B_k_plus;
  emxArray_real_T *Delta;
  emxArray_real_T *S_k;
  emxArray_real_T *d_k;
  real_T(*x_ref)[105];
  real_T(*u_ref)[30];
  real_T(*tspan)[2];
  real_T L;
  real_T N;
  real_T N_sub;
  real_T alpha;
  real_T b_I;
  real_T s_ref;
  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  /* Marshall function inputs */
  N = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "N");
  tspan = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "tspan");
  x_ref = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "x_ref");
  u_ref = g_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "u_ref");
  s_ref = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "s_ref");
  N_sub = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "N_sub");
  L = emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "L");
  b_I = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "I");
  alpha = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "alpha");
  /* Invoke the target function */
  emxInit_real_T(&st, &A_k, 3, &id_emlrtRTEI);
  emxInit_real_T(&st, &B_k_plus, 3, &id_emlrtRTEI);
  emxInit_real_T(&st, &B_k_minus, 3, &id_emlrtRTEI);
  emxInit_real_T(&st, &S_k, 3, &id_emlrtRTEI);
  emxInit_real_T(&st, &d_k, 3, &id_emlrtRTEI);
  emxInit_real_T(&st, &Delta, 2, &id_emlrtRTEI);
  discretize_error_dynamics_FOH_RKV65_3DoF(
      &st, N, *tspan, *x_ref, *u_ref, s_ref, N_sub, L, b_I, alpha, A_k,
      B_k_plus, B_k_minus, S_k, d_k, Delta);
  /* Marshall function outputs */
  A_k->canFreeData = false;
  plhs[0] = emlrt_marshallOut(A_k);
  emxFree_real_T(&st, &A_k);
  if (nlhs > 1) {
    B_k_plus->canFreeData = false;
    plhs[1] = emlrt_marshallOut(B_k_plus);
  }
  emxFree_real_T(&st, &B_k_plus);
  if (nlhs > 2) {
    B_k_minus->canFreeData = false;
    plhs[2] = emlrt_marshallOut(B_k_minus);
  }
  emxFree_real_T(&st, &B_k_minus);
  if (nlhs > 3) {
    S_k->canFreeData = false;
    plhs[3] = emlrt_marshallOut(S_k);
  }
  emxFree_real_T(&st, &S_k);
  if (nlhs > 4) {
    d_k->canFreeData = false;
    plhs[4] = emlrt_marshallOut(d_k);
  }
  emxFree_real_T(&st, &d_k);
  if (nlhs > 5) {
    Delta->canFreeData = false;
    plhs[5] = b_emlrt_marshallOut(Delta);
  }
  emxFree_real_T(&st, &Delta);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation
 * (_coder_discretize_error_dynamics_FOH_RKV65_3DoF_api.c) */
