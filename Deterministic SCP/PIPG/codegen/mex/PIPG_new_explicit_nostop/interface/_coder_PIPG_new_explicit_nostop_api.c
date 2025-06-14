/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PIPG_new_explicit_nostop_api.c
 *
 * Code generation for function '_coder_PIPG_new_explicit_nostop_api'
 *
 */

/* Include files */
#include "_coder_PIPG_new_explicit_nostop_api.h"
#include "PIPG_new_explicit_nostop.h"
#include "PIPG_new_explicit_nostop_data.h"
#include "PIPG_new_explicit_nostop_emxutil.h"
#include "PIPG_new_explicit_nostop_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo g_emlrtRTEI = {
    1,                                     /* lineNo */
    1,                                     /* colNo */
    "_coder_PIPG_new_explicit_nostop_api", /* fName */
    ""                                     /* pName */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[1000];

static const mxArray *b_emlrt_marshallOut(real_T u[3]);

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3000];

static const mxArray *c_emlrt_marshallOut(const emlrtStack *sp,
                                          const struct0_T u);

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3000];

static const mxArray *d_emlrt_marshallOut(const emxArray_real_T *u);

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3];

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier))[1000];

static const mxArray *emlrt_marshallOut(real_T u[1000]);

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3];

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[1000000];

static real_T (
    *h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                        const emlrtMsgIdentifier *parentId))[1000000];

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier);

static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[1000];

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3000];

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3];

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[1000000];

static real_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[1000]
{
  real_T(*y)[1000];
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *b_emlrt_marshallOut(real_T u[3])
{
  static const int32_T i = 0;
  static const int32_T i1 = 3;
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u[0]);
  emlrtSetDimensions((mxArray *)m, &i1, 1);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[3000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *c_emlrt_marshallOut(const emlrtStack *sp,
                                          const struct0_T u)
{
  static const int32_T iv[2] = {1, 8};
  static const char_T *sv[8] = {"primal_infeasible",
                                "dual_infeasible",
                                "solution_status",
                                "iterations",
                                "zhat_inf_dj",
                                "what_inf_dj",
                                "zhis",
                                "time"};
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *m;
  const mxArray *propValues;
  const mxArray *y;
  real_T *pData;
  int32_T iv1[3];
  int32_T b_i;
  int32_T c_i;
  int32_T d_i;
  int32_T i;
  const char_T *propClasses = "coder.internal.string";
  const char_T *propNames = "Value";
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 8, (const char_T **)&sv[0]));
  b_y = NULL;
  propValues = emlrtCreateLogicalScalar(u.primal_infeasible);
  emlrtAssign(&b_y, propValues);
  emlrtSetFieldR2017b(y, 0, "primal_infeasible", b_y, 0);
  c_y = NULL;
  propValues = emlrtCreateLogicalScalar(u.dual_infeasible);
  emlrtAssign(&c_y, propValues);
  emlrtSetFieldR2017b(y, 0, "dual_infeasible", c_y, 1);
  d_y = NULL;
  m = NULL;
  emlrtAssign(&d_y, emlrtCreateClassInstance2022a((emlrtCTX)sp,
                                                  "coder.internal.string"));
  m = NULL;
  e_y = NULL;
  propValues = emlrtCreateCharArray(2, &iv[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 8, propValues,
                           &u.solution_status.Value[0]);
  emlrtAssign(&e_y, propValues);
  emlrtAssign(&m, e_y);
  propValues = m;
  emlrtSetAllProperties((emlrtCTX)sp, &d_y, 0, 1, (const char_T **)&propNames,
                        (const char_T **)&propClasses, &propValues);
  emlrtAssign(&d_y, emlrtConvertInstanceToRedirectSource(
                        (emlrtCTX)sp, d_y, 0, "coder.internal.string"));
  emlrtSetFieldR2017b(y, 0, "solution_status", d_y, 2);
  f_y = NULL;
  propValues = emlrtCreateDoubleScalar(u.iterations);
  emlrtAssign(&f_y, propValues);
  emlrtSetFieldR2017b(y, 0, "iterations", f_y, 3);
  emlrtSetFieldR2017b(y, 0, "zhat_inf_dj", d_emlrt_marshallOut(u.zhat_inf_dj),
                      4);
  emlrtSetFieldR2017b(y, 0, "what_inf_dj", d_emlrt_marshallOut(u.what_inf_dj),
                      5);
  g_y = NULL;
  iv1[0] = u.zhis->size[0];
  iv1[1] = u.zhis->size[1];
  iv1[2] = u.zhis->size[2];
  propValues = emlrtCreateNumericArray(3, &iv1[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(propValues);
  i = 0;
  for (b_i = 0; b_i < u.zhis->size[2]; b_i++) {
    for (c_i = 0; c_i < u.zhis->size[1]; c_i++) {
      for (d_i = 0; d_i < u.zhis->size[0]; d_i++) {
        pData[i] = u.zhis->data[(d_i + u.zhis->size[0] * c_i) +
                                u.zhis->size[0] * u.zhis->size[1] * b_i];
        i++;
      }
    }
  }
  emlrtAssign(&g_y, propValues);
  emlrtSetFieldR2017b(y, 0, "zhis", g_y, 6);
  h_y = NULL;
  propValues = emlrtCreateDoubleScalar(u.time);
  emlrtAssign(&h_y, propValues);
  emlrtSetFieldR2017b(y, 0, "time", h_y, 7);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3000]
{
  real_T(*y)[3000];
  y = l_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *d_emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *m;
  const mxArray *y;
  const real_T *u_data;
  real_T *pData;
  int32_T b_i;
  int32_T i;
  u_data = u->data;
  y = NULL;
  m = emlrtCreateNumericArray(1, &u->size[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < u->size[0]; b_i++) {
    pData[i] = u_data[b_i];
    i++;
  }
  emlrtAssign(&y, m);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[3];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier))[1000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[1000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *emlrt_marshallOut(real_T u[1000])
{
  static const int32_T i = 0;
  static const int32_T i1 = 1000;
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u[0]);
  emlrtSetDimensions((mxArray *)m, &i1, 1);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3]
{
  real_T(*y)[3];
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[1000000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[1000000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[1000000]
{
  real_T(*y)[1000000];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[1000]
{
  static const int32_T dims = 1000;
  real_T(*ret)[1000];
  int32_T i;
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real_T(*)[1000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3000]
{
  static const int32_T dims[2] = {3, 1000};
  real_T(*ret)[3000];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[3000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3]
{
  static const int32_T dims = 3;
  real_T(*ret)[3];
  int32_T i;
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real_T(*)[3])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[1000000]
{
  static const int32_T dims[2] = {1000, 1000};
  real_T(*ret)[1000000];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[1000000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

void PIPG_new_explicit_nostop_api(const mxArray *const prhs[16], int32_T nlhs,
                                  const mxArray *plhs[3])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  struct0_T sol_info;
  real_T(*L)[1000000];
  real_T(*L_inv)[1000000];
  real_T(*Hhat)[3000];
  real_T(*qhat)[1000];
  real_T(*z_ref)[1000];
  real_T(*z_star)[1000];
  real_T(*h)[3];
  real_T(*what_ref)[3];
  real_T(*what_star)[3];
  real_T j_check;
  real_T j_max;
  real_T lambda;
  real_T omega;
  real_T rho;
  real_T sigma;
  real_T tol_abs;
  real_T tol_infeas;
  real_T tol_rel;
  st.tls = emlrtRootTLSGlobal;
  z_star = (real_T(*)[1000])mxMalloc(sizeof(real_T[1000]));
  what_star = (real_T(*)[3])mxMalloc(sizeof(real_T[3]));
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  /* Marshall function inputs */
  qhat = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "qhat");
  Hhat = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "Hhat");
  h = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "h");
  L = g_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "L");
  L_inv = g_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "L_inv");
  lambda = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "lambda");
  sigma = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "sigma");
  omega = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "omega");
  rho = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "rho");
  tol_abs = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "tol_abs");
  tol_rel = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "tol_rel");
  tol_infeas = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "tol_infeas");
  j_check = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[12]), "j_check");
  j_max = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[13]), "j_max");
  z_ref = emlrt_marshallIn(&st, emlrtAlias(prhs[14]), "z_ref");
  what_ref = e_emlrt_marshallIn(&st, emlrtAlias(prhs[15]), "what_ref");
  /* Invoke the target function */
  emxInitStruct_struct0_T(&st, &sol_info, &g_emlrtRTEI);
  PIPG_new_explicit_nostop(&st, *qhat, *Hhat, *h, *L, *L_inv, lambda, sigma,
                           omega, rho, tol_abs, tol_rel, tol_infeas, j_check,
                           j_max, *z_ref, *what_ref, *z_star, *what_star,
                           &sol_info);
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*z_star);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(*what_star);
  }
  if (nlhs > 2) {
    plhs[2] = c_emlrt_marshallOut(&st, sol_info);
  }
  emxFreeStruct_struct0_T(&st, &sol_info);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_PIPG_new_explicit_nostop_api.c) */
