/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PIPG_orig_precond_Dexplicit_api.c
 *
 * Code generation for function '_coder_PIPG_orig_precond_Dexplicit_api'
 *
 */

/* Include files */
#include "_coder_PIPG_orig_precond_Dexplicit_api.h"
#include "PIPG_orig_precond_Dexplicit.h"
#include "PIPG_orig_precond_Dexplicit_data.h"
#include "PIPG_orig_precond_Dexplicit_emxutil.h"
#include "PIPG_orig_precond_Dexplicit_mexutil.h"
#include "PIPG_orig_precond_Dexplicit_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo r_emlrtRTEI = {
    1,                                        /* lineNo */
    1,                                        /* colNo */
    "_coder_PIPG_orig_precond_Dexplicit_api", /* fName */
    ""                                        /* pName */
};

/* Function Declarations */
static int32_T ab_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId);

static const mxArray *b_emlrt_marshallOut(real_T u[3]);

static real_T (*bb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId))[3];

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[1000];

static const mxArray *c_emlrt_marshallOut(const emlrtStack *sp,
                                          const struct0_T u);

static void cb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[2]);

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[1000];

static const mxArray *d_emlrt_marshallOut(const emxArray_real_T *u);

static void db_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId);

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, sparse *y);

static void eb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[3]);

static const mxArray *emlrt_marshallOut(real_T u[1000]);

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, sparse *y);

static void fb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[3]);

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static real_T (*gb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId))[1000000];

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_int32_T *y);

static int32_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                  const emlrtMsgIdentifier *parentId);

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3];

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3];

static Ballz l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                const char_T *identifier);

static Ballz m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                const emlrtMsgIdentifier *parentId);

static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[2]);

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId);

static Singletonz p_emlrt_marshallIn(const emlrtStack *sp,
                                     const mxArray *nullptr,
                                     const char_T *identifier);

static Singletonz q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                     const emlrtMsgIdentifier *parentId);

static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3]);

static void s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3]);

static real_T (*t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[1000000];

static real_T (
    *u_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                        const emlrtMsgIdentifier *parentId))[1000000];

static real_T (*w_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[1000];

static void x_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void y_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_int32_T *ret);

/* Function Definitions */
static int32_T ab_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  int32_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "int32", false, 0U,
                          (const void *)&dims);
  ret = *(int32_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
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
  static const char_T *sv[7] = {"primal_infeasible",
                                "dual_infeasible",
                                "solution_status",
                                "iterations",
                                "zhat_inf_dj",
                                "what_inf_dj",
                                "zhis"};
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *g_y;
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
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 7, (const char_T **)&sv[0]));
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
  return y;
}

static void cb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[1000]
{
  real_T(*y)[1000];
  y = w_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
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

static void db_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims[2] = {0, 0};
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                          (const void *)&dims[0]);
  emlrtMxGetData(src);
  emlrtDestroyArray(&src);
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, sparse *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId, y);
  emlrtDestroyArray(&nullptr);
}

static void eb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, sparse *y)
{
  emlrtMsgIdentifier thisId;
  const mxArray *propValues[4];
  int32_T iv[2];
  const char_T *propClasses[4] = {
      "coder.internal.sparse", "coder.internal.sparse", "coder.internal.sparse",
      "coder.internal.sparse"};
  const char_T *propNames[4] = {"d", "colidx", "rowidx", "maxnz"};
  boolean_T bv[2];
  propValues[0] = NULL;
  propValues[1] = NULL;
  propValues[2] = NULL;
  propValues[3] = NULL;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  bv[0] = false;
  iv[0] = 3;
  bv[1] = false;
  iv[1] = 1000;
  emlrtCheckSparse((emlrtCTX)sp, parentId, u, &iv[0], &bv[0], mxDOUBLE_CLASS,
                   mxREAL);
  emlrtCheckMcosClass2017a((emlrtCTX)sp, parentId, u, "sparse");
  emlrtAssign(&u, emlrtConvertInstanceToRedirectTarget(
                      (emlrtCTX)sp, u, 0, "coder.internal.sparse"));
  emlrtCheckMcosClass2017a((emlrtCTX)sp, parentId, u, "coder.internal.sparse");
  emlrtGetAllProperties((emlrtCTX)sp, u, 0, 4, (const char_T **)&propNames[0],
                        (const char_T **)&propClasses[0], &propValues[0]);
  thisId.fIdentifier = "d";
  g_emlrt_marshallIn(sp, emlrtAlias(propValues[0]), &thisId, y->d);
  thisId.fIdentifier = "colidx";
  h_emlrt_marshallIn(sp, emlrtAlias(propValues[1]), &thisId, y->colidx);
  thisId.fIdentifier = "rowidx";
  h_emlrt_marshallIn(sp, emlrtAlias(propValues[2]), &thisId, y->rowidx);
  thisId.fIdentifier = "maxnz";
  i_emlrt_marshallIn(sp, emlrtAlias(propValues[3]), &thisId);
  emlrtDestroyArrays(4, &propValues[0]);
  emlrtDestroyArray(&u);
}

static void fb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  x_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T (*gb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_int32_T *y)
{
  y_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static int32_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                  const emlrtMsgIdentifier *parentId)
{
  int32_T y;
  y = ab_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[3];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = k_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3]
{
  real_T(*y)[3];
  y = bb_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static Ballz l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                const char_T *identifier)
{
  Ballz y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = m_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static Ballz m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
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
  n_emlrt_marshallIn(sp, emlrtAlias(propValues[0]), &thisId, y.indices);
  thisId.fIdentifier = "k_override";
  o_emlrt_marshallIn(sp, emlrtAlias(propValues[1]), &thisId);
  thisId.fIdentifier = "r";
  y.r = b_emlrt_marshallIn(sp, emlrtAlias(propValues[2]), &thisId);
  emlrtDestroyArrays(3, &propValues[0]);
  emlrtDestroyArray(&u);
  return y;
}

static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[2])
{
  cb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId)
{
  db_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
}

static Singletonz p_emlrt_marshallIn(const emlrtStack *sp,
                                     const mxArray *nullptr,
                                     const char_T *identifier)
{
  Singletonz y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = q_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static Singletonz q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
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
  r_emlrt_marshallIn(sp, emlrtAlias(propValues[0]), &thisId, y.indices);
  thisId.fIdentifier = "k_override";
  o_emlrt_marshallIn(sp, emlrtAlias(propValues[1]), &thisId);
  thisId.fIdentifier = "y";
  s_emlrt_marshallIn(sp, emlrtAlias(propValues[2]), &thisId, y.y);
  emlrtDestroyArrays(3, &propValues[0]);
  emlrtDestroyArray(&u);
  return y;
}

static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3])
{
  eb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3])
{
  fb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T (*t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[1000000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[1000000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = u_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*u_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[1000000]
{
  real_T(*y)[1000000];
  y = gb_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*w_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static void x_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims = -1;
  real_T *ret_data;
  int32_T i;
  int32_T i1;
  boolean_T b = true;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  i1 = ret->size[0];
  ret->size[0] = i;
  emxEnsureCapacity_real_T(sp, ret, i1, (emlrtRTEInfo *)NULL);
  ret_data = ret->data;
  emlrtImportArrayR2015b((emlrtConstCTX)sp, src, &ret_data[0], 8, false);
  emlrtDestroyArray(&src);
}

static void y_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_int32_T *ret)
{
  static const int32_T dims = -1;
  int32_T i;
  int32_T i1;
  int32_T *ret_data;
  boolean_T b = true;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "int32", false, 1U,
                            (const void *)&dims, &b, &i);
  i1 = ret->size[0];
  ret->size[0] = i;
  emxEnsureCapacity_int32_T(sp, ret, i1, (emlrtRTEInfo *)NULL);
  ret_data = ret->data;
  emlrtImportArrayR2015b((emlrtConstCTX)sp, src, &ret_data[0], 4, false);
  emlrtDestroyArray(&src);
}

void PIPG_orig_precond_Dexplicit_api(c_PIPG_orig_precond_DexplicitSt *SD,
                                     const mxArray *const prhs[20],
                                     int32_T nlhs, const mxArray *plhs[3])
{
  Ballz Dball;
  Singletonz Dsingleton;
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  sparse Hhat;
  struct0_T sol_info;
  real_T(*L)[1000000];
  real_T(*L_inv)[1000000];
  real_T(*qhat)[1000];
  real_T(*z_ref)[1000];
  real_T(*z_star)[1000];
  real_T(*h)[3];
  real_T(*what_ref)[3];
  real_T(*what_star)[3];
  real_T S_z;
  real_T c_z;
  real_T j_check;
  real_T j_max;
  real_T j_restart;
  real_T lambda;
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
  emxInitStruct_sparse(&st, &Hhat, &r_emlrtRTEI);
  e_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "Hhat", &Hhat);
  h = j_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "h");
  Dball = l_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "Dball");
  Dsingleton = p_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "Dsingleton");
  L = t_emlrt_marshallIn(&st, emlrtAlias(prhs[5]), "L");
  L_inv = t_emlrt_marshallIn(&st, emlrtAlias(prhs[6]), "L_inv");
  lambda = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "lambda");
  sigma = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "sigma");
  rho = emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "rho");
  tol_abs = emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "tol_abs");
  tol_rel = emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "tol_rel");
  tol_infeas = emlrt_marshallIn(&st, emlrtAliasP(prhs[12]), "tol_infeas");
  j_check = emlrt_marshallIn(&st, emlrtAliasP(prhs[13]), "j_check");
  j_restart = emlrt_marshallIn(&st, emlrtAliasP(prhs[14]), "j_restart");
  j_max = emlrt_marshallIn(&st, emlrtAliasP(prhs[15]), "j_max");
  z_ref = c_emlrt_marshallIn(&st, emlrtAlias(prhs[16]), "z_ref");
  what_ref = j_emlrt_marshallIn(&st, emlrtAlias(prhs[17]), "what_ref");
  S_z = emlrt_marshallIn(&st, emlrtAliasP(prhs[18]), "S_z");
  c_z = emlrt_marshallIn(&st, emlrtAliasP(prhs[19]), "c_z");
  /* Invoke the target function */
  emxInitStruct_struct0_T(&st, &sol_info, &r_emlrtRTEI);
  PIPG_orig_precond_Dexplicit(
      SD, &st, *qhat, Hhat, *h, &Dball, &Dsingleton, *L, *L_inv, lambda, sigma,
      rho, tol_abs, tol_rel, tol_infeas, j_check, j_restart, j_max, *z_ref,
      *what_ref, S_z, c_z, *z_star, *what_star, &sol_info);
  emxFreeStruct_sparse(&st, &Hhat);
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

/* End of code generation (_coder_PIPG_orig_precond_Dexplicit_api.c) */
