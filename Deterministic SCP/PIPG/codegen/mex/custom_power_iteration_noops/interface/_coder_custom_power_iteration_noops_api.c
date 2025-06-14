/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_custom_power_iteration_noops_api.c
 *
 * Code generation for function '_coder_custom_power_iteration_noops_api'
 *
 */

/* Include files */
#include "_coder_custom_power_iteration_noops_api.h"
#include "custom_power_iteration_noops.h"
#include "custom_power_iteration_noops_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[686];

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[196];

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[196];

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[98];

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier))[686];

static const mxArray *emlrt_marshallOut(const real_T u);

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[98];

static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier);

static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[105];

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[105];

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[30];

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[30];

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[686];

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[196];

static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[98];

static real_T p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[105];

static real_T (*r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[30];

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[686]
{
  real_T(*y)[686];
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[196]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[196];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[196]
{
  real_T(*y)[196];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[98]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[98];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier))[686]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[686];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[98]
{
  real_T(*y)[98];
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = p_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[105]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[105];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[105]
{
  real_T(*y)[105];
  y = q_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[30]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[30];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = l_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[30]
{
  real_T(*y)[30];
  y = r_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[686]
{
  static const int32_T dims[3] = {7, 7, 14};
  real_T(*ret)[686];
  int32_T iv[3];
  boolean_T bv[3] = {false, false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 3U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[686])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[196]
{
  static const int32_T dims[3] = {7, 2, 14};
  real_T(*ret)[196];
  int32_T iv[3];
  boolean_T bv[3] = {false, false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 3U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[196])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[98]
{
  static const int32_T dims[3] = {7, 1, 14};
  real_T(*ret)[98];
  int32_T iv[3];
  boolean_T bv[3] = {false, false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 3U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[98])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static real_T (*r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

void e_custom_power_iteration_noops_(const mxArray *const prhs[11],
                                     const mxArray **plhs)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *prhs_copy_idx_7;
  const mxArray *prhs_copy_idx_8;
  const mxArray *prhs_copy_idx_9;
  real_T(*Ahat_1)[686];
  real_T(*Ahat_2)[686];
  real_T(*Bhat_minus)[196];
  real_T(*Bhat_plus)[196];
  real_T(*x)[105];
  real_T(*xi)[105];
  real_T(*Shat)[98];
  real_T(*u)[30];
  real_T lx1_inv;
  real_T lx2_inv;
  real_T s;
  st.tls = emlrtRootTLSGlobal;
  prhs_copy_idx_7 = emlrtProtectR2012b(prhs[7], 7, false, -1);
  prhs_copy_idx_8 = emlrtProtectR2012b(prhs[8], 8, false, -1);
  prhs_copy_idx_9 = emlrtProtectR2012b(prhs[9], 9, false, -1);
  /* Marshall function inputs */
  Ahat_1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "Ahat_1");
  Ahat_2 = emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "Ahat_2");
  Bhat_minus = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "Bhat_minus");
  Bhat_plus = c_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "Bhat_plus");
  Shat = e_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "Shat");
  lx1_inv = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "lx1_inv");
  lx2_inv = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "lx2_inv");
  x = i_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_7), "x");
  xi = i_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_8), "xi");
  u = k_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_9), "u");
  s = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "s");
  /* Invoke the target function */
  lx1_inv = custom_power_iteration_noops(&st, *Ahat_1, *Ahat_2, *Bhat_minus,
                                         *Bhat_plus, *Shat, lx1_inv, lx2_inv,
                                         *x, *xi, *u, s);
  /* Marshall function outputs */
  *plhs = emlrt_marshallOut(lx1_inv);
}

/* End of code generation (_coder_custom_power_iteration_noops_api.c) */
