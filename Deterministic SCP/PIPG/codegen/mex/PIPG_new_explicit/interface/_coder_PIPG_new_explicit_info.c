/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PIPG_new_explicit_info.c
 *
 * Code generation for function 'PIPG_new_explicit'
 *
 */

/* Include files */
#include "_coder_PIPG_new_explicit_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[7] = {
      "789ced984d8fd24018c707b3ba1b139583eec9035f808dec8bc97a630b0b9bb8d0f0b21a"
      "a8596a3b2c85b6039de1f50be8cde845132fde3dfb094cfc0826c69b"
      "7bf123ecd1423b1426694a000bcbf21c983efca7cf6fe6a1f933000227a70100c07d60c5"
      "a7c7d678cfce83f6780b8c07ab07ecf1369383e1fb1b63f751fd9d3d",
      "4a4827b043ac44173538bc53469aa28b3ac975eb10181023b505e581525654985334981d"
      "4d52fd4c3b1e9186495fea5f731528d5b24d0d1815ecac501d4d86fd"
      "28b9ec77c3a31f6cb0fd60e7cdcaf3ea7fd0568e4455edcd81c7c6a4fbbb9a9247eb7ff0"
      "e051bdc8bde29e09790c0d2c908a48ca425a87314369412186a4a606",
      "758205be69c84d18e2541163888568341e3a387c12e5a2fd4beb2a9c2548aa8898285238"
      "cbf1e10c926a90849f8bbaace817420c126898cfa7d29f10322708fc"
      "099f1038a4b76027640e66d941cb77b4f13e945cf6f960c23eb0a3337f6b307efdf93be0"
      "27ef217e73e9278fc6a2781d977a933ec7db2ebc20a3772238b29fde",
      "6d4773472f729dbd33dcab2991b8b30ede83e3b50ee092fb55bfe472fffff29ffa943c5a"
      "ffae078fea59d31d544890de1be797a6e4b3e197df7ef4e0517d89fc"
      "d669bd69baabeab79bdf76affce4d15875bf2da7cdb3e74b24658e51b2dda816aa8d4c24"
      "c1adfd960dbfceb39bae3c4ba91ba80a25b2b0f3339d382defb32b6f",
      "5c2fa697c760ed036d80f61e809be343abcef3cb67f385d4d9d34436564d251bf9ea4597"
      "17f77829b1f6593696d567177d8e9ed577bf78f0a9be4cbe3b7ab0a5"
      "e67b537c69d5796bdf9d4ffd597fefbe77a91f64f4627c91b6d07f39d761fb1c76eaaa22"
      "2984febd38b7ef832d2677fa602998a07add5ce1bc78775c799622a3",
      "e66b153afbfb3b23efad2b6f5c5fece74cbbbca3f9e54bbf7ef8fb3fc6e5f73fa77ef268"
      "5c579f7fe4c20b32ba46a4d8c121ce2722fb8554a65b8918491983eb"
      "eff3ff00d653e157",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 7344U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *propFieldName[9] = {"Version",
                                    "ResolvedFunctions",
                                    "Checksum",
                                    "EntryPoints",
                                    "CoverageInfo",
                                    "IsPolymorphic",
                                    "PropertyList",
                                    "UUID",
                                    "ClassEntryPointIsHandle"};
  const char_T *epFieldName[8] = {
      "QualifiedName",    "NumberOfInputs", "NumberOfOutputs", "ConstantInputs",
      "ResolvedFilePath", "TimeStamp",      "Constructor",     "Visible"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 1, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 18);
  emlrtSetField(xEntryPoints, 0, "QualifiedName",
                emlrtMxCreateString("PIPG_new_explicit"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(18.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
          "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic "
          "SCP\\P"
          "IPG\\PIPG_new_explicit.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739775.63859953708));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("24.2.0.2740171 (R2024b) Update 1"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("Vq0KlPQNhgSAgsUrrYlkIB"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_PIPG_new_explicit_info.c) */
