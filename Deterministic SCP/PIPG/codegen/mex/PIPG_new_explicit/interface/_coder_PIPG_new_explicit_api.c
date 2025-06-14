/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PIPG_new_explicit_api.c
 *
 * Code generation for function '_coder_PIPG_new_explicit_api'
 *
 */

/* Include files */
#include "_coder_PIPG_new_explicit_api.h"
#include "PIPG_new_explicit.h"
#include "PIPG_new_explicit_data.h"
#include "PIPG_new_explicit_emxutil.h"
#include "PIPG_new_explicit_mexutil.h"
#include "PIPG_new_explicit_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo h_emlrtRTEI = {
    1,                              /* lineNo */
    1,                              /* colNo */
    "_coder_PIPG_new_explicit_api", /* fName */
    ""                              /* pName */
};

/* Function Declarations */
static void ab_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[3]);

static const mxArray *b_emlrt_marshallOut(real_T u[3]);

static real_T (*bb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId))[1000000];

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[1000];

static const mxArray *c_emlrt_marshallOut(const emlrtStack *sp,
                                          const struct0_T u);

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[1000];

static const mxArray *d_emlrt_marshallOut(const emxArray_real_T *u);

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3000];

static const mxArray *emlrt_marshallOut(real_T u[1000]);

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3000];

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3];

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3];

static Ballz i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                const char_T *identifier);

static Ballz j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                const emlrtMsgIdentifier *parentId);

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[2]);

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId);

static Singletonz m_emlrt_marshallIn(const emlrtStack *sp,
                                     const mxArray *nullptr,
                                     const char_T *identifier);

static Singletonz n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                     const emlrtMsgIdentifier *parentId);

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3]);

static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3]);

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[1000000];

static real_T (
    *r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                        const emlrtMsgIdentifier *parentId))[1000000];

static real_T (*t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[1000];

static real_T (*u_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3000];

static real_T (*v_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3];

static void w_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[2]);

static void x_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId);

static void y_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[3]);

/* Function Definitions */
static void ab_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[3])
{
  static const int32_T dims = 3;
  real_T(*r)[3];
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                          (const void *)&dims);
  r = (real_T(*)[3])emlrtMxGetData(src);
  ret[0] = (*r)[0];
  ret[1] = (*r)[1];
  ret[2] = (*r)[2];
  emlrtDestroyArray(&src);
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

static real_T (*bb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[1000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[1000];
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
  int32_T iv[2];
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
  iv[0] = 1;
  iv[1] = u.solution_status.Value.size[1];
  propValues = emlrtCreateCharArray(2, &iv[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, u.solution_status.Value.size[1],
                           propValues, &u.solution_status.Value.data[0]);
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
                                   const emlrtMsgIdentifier *parentId))[1000]
{
  real_T(*y)[1000];
  y = t_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
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
                                   const char_T *identifier))[3000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[3000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
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
                                   const emlrtMsgIdentifier *parentId))[3000]
{
  real_T(*y)[3000];
  y = u_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[3];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3]
{
  real_T(*y)[3];
  y = v_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static Ballz i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                const char_T *identifier)
{
  Ballz y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static Ballz j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                const emlrtMsgIdentifier *parentId)
{
  Ballz y;
  emlrtMsgIdentifier thisId;
  const mxArray *propValues[3];
  const char_T *propClasses[3] = {"Ballz", "Ballz", "Ballz"};
  const char_T *propNames[3] = {"indices", "k_override", "r"};
  propValues[0] = NULL;
  propValues[1] = NULL;
  propValues[2] = NULL;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckMcosClass2017a((emlrtCTX)sp, parentId, u, "Ballz");
  emlrtGetAllProperties((emlrtCTX)sp, u, 0, 3, (const char_T **)&propNames[0],
                        (const char_T **)&propClasses[0], &propValues[0]);
  thisId.fIdentifier = "indices";
  k_emlrt_marshallIn(sp, emlrtAlias(propValues[0]), &thisId, y.indices);
  thisId.fIdentifier = "k_override";
  l_emlrt_marshallIn(sp, emlrtAlias(propValues[1]), &thisId);
  thisId.fIdentifier = "r";
  y.r = b_emlrt_marshallIn(sp, emlrtAlias(propValues[2]), &thisId);
  emlrtDestroyArrays(3, &propValues[0]);
  emlrtDestroyArray(&u);
  return y;
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[2])
{
  w_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId)
{
  x_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
}

static Singletonz m_emlrt_marshallIn(const emlrtStack *sp,
                                     const mxArray *nullptr,
                                     const char_T *identifier)
{
  Singletonz y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = n_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static Singletonz n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                     const emlrtMsgIdentifier *parentId)
{
  Singletonz y;
  emlrtMsgIdentifier thisId;
  const mxArray *propValues[3];
  const char_T *propClasses[3] = {"Singletonz", "Singletonz", "Singletonz"};
  const char_T *propNames[3] = {"indices", "k_override", "y"};
  propValues[0] = NULL;
  propValues[1] = NULL;
  propValues[2] = NULL;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckMcosClass2017a((emlrtCTX)sp, parentId, u, "Singletonz");
  emlrtGetAllProperties((emlrtCTX)sp, u, 0, 3, (const char_T **)&propNames[0],
                        (const char_T **)&propClasses[0], &propValues[0]);
  thisId.fIdentifier = "indices";
  o_emlrt_marshallIn(sp, emlrtAlias(propValues[0]), &thisId, y.indices);
  thisId.fIdentifier = "k_override";
  l_emlrt_marshallIn(sp, emlrtAlias(propValues[1]), &thisId);
  thisId.fIdentifier = "y";
  p_emlrt_marshallIn(sp, emlrtAlias(propValues[2]), &thisId, y.y);
  emlrtDestroyArrays(3, &propValues[0]);
  emlrtDestroyArray(&u);
  return y;
}

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3])
{
  y_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3])
{
  ab_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[1000000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[1000000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = r_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[1000000]
{
  real_T(*y)[1000000];
  y = bb_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static real_T (*u_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static real_T (*v_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static void w_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[2])
{
  static const int32_T dims[2] = {1, 2};
  real_T(*r)[2];
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                          (const void *)&dims[0]);
  r = (real_T(*)[2])emlrtMxGetData(src);
  ret[0] = (*r)[0];
  ret[1] = (*r)[1];
  emlrtDestroyArray(&src);
}

static void x_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims[2] = {0, 0};
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                          (const void *)&dims[0]);
  emlrtMxGetData(src);
  emlrtDestroyArray(&src);
}

static void y_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[3])
{
  static const int32_T dims[2] = {1, 3};
  real_T(*r)[3];
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                          (const void *)&dims[0]);
  r = (real_T(*)[3])emlrtMxGetData(src);
  ret[0] = (*r)[0];
  ret[1] = (*r)[1];
  ret[2] = (*r)[2];
  emlrtDestroyArray(&src);
}

void PIPG_new_explicit_api(const mxArray *const prhs[18], int32_T nlhs,
                           const mxArray *plhs[3])
{
  Ballz Dball;
  Singletonz Dsingleton;
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
  qhat = c_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "qhat");
  Hhat = e_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "Hhat");
  h = g_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "h");
  Dball = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "Dball");
  Dsingleton = m_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "Dsingleton");
  L = q_emlrt_marshallIn(&st, emlrtAlias(prhs[5]), "L");
  L_inv = q_emlrt_marshallIn(&st, emlrtAlias(prhs[6]), "L_inv");
  lambda = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "lambda");
  sigma = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "sigma");
  omega = emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "omega");
  rho = emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "rho");
  tol_abs = emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "tol_abs");
  tol_rel = emlrt_marshallIn(&st, emlrtAliasP(prhs[12]), "tol_rel");
  tol_infeas = emlrt_marshallIn(&st, emlrtAliasP(prhs[13]), "tol_infeas");
  j_check = emlrt_marshallIn(&st, emlrtAliasP(prhs[14]), "j_check");
  j_max = emlrt_marshallIn(&st, emlrtAliasP(prhs[15]), "j_max");
  z_ref = c_emlrt_marshallIn(&st, emlrtAlias(prhs[16]), "z_ref");
  what_ref = g_emlrt_marshallIn(&st, emlrtAlias(prhs[17]), "what_ref");
  /* Invoke the target function */
  emxInitStruct_struct0_T(&st, &sol_info, &h_emlrtRTEI);
  PIPG_new_explicit(&st, *qhat, *Hhat, *h, &Dball, &Dsingleton, *L, *L_inv,
                    lambda, sigma, omega, rho, tol_abs, tol_rel, tol_infeas,
                    j_check, j_max, *z_ref, *what_ref, *z_star, *what_star,
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

/* End of code generation (_coder_PIPG_new_explicit_api.c) */
