/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_discretize_error_dynamics_FOH_RKV65_3DoF_info.c
 *
 * Code generation for function 'discretize_error_dynamics_FOH_RKV65_3DoF'
 *
 */

/* Include files */
#include "_coder_discretize_error_dynamics_FOH_RKV65_3DoF_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[8] = {
      "789ced993d6fd3401cc62fa820182819104248886c744955681b5e06a4344e52ab4993c6"
      "498bc0c875ed4b7dadcfe7ded96d92890131b074841dd8d9f80888a5"
      "1b2c7c01464646e226ce8b85d548491d487ccbf99fe7ece7fcdce97eb202227c3e020098"
      "05adf6e94eabbfdaaea3edfe02e86f5e3dd2ee2f7a6ad0f97da6ef3e",
      "573f6ef70a312c58b35a852163d8b953251819b26195eb26041432a21f42f554a9221d96"
      "1186426fb1ee5438d323750a4772ae531a54f6051b03aab1ee0cf5de"
      "a293c75ce4efef3b33601e5f7cf2887af4e7e917a9c7628541ca444b93adaa58302047d1"
      "211439a2d8181a16138b36556d184be932639089c9643ab6fc682199",
      "4a3a97adabb860114593998594b8902ac64b44d987563c271b2a3276450e5a9036f344ce"
      "8058738058e48b5991434ca1d0420dd942c41055b78412a4945049ad"
      "379704294cca1456a5d2da6662595ae448661e4784725e5251b52ac10347ece4b6ed93cb"
      "a0b95df2cdada5244f270002f35389bda3c3aedfcf21fd5ef9faf5eb",
      "41ed8bf602c7f244853a139d70c5647b919d799c95f3b501dfdbdb77c75f3eed3f7cff11"
      "09d2af7072fd6d907e6e1b975fcde77983eedb1b3e7e518fae6e2435"
      "6ebd71ff68ef61dae4c9f2bdacb6c4af76e7513cc3e7ac79009f3aa8e7875cf837b9b012"
      "72e1dcb9b032055c787df776c805307a2ea0fce262a22eb30a9f4869",
      "f93d212df0b9dd74c88569e7c24b9f5c06cdedd619b9b9ba50c7ee99e6cc48c2cd682493"
      "e8326de9e3e2c6ef21fddef9faf5ebe3e4867ff6f378523932fbf54a"
      "c811307a8e60ad56e1d1026fd20dc4d5682a71806c940d3932ed1cd9f6c96554e7b6107e"
      "5f9c3f27a6e0fb62e7fd5cc805307a2e6ce5cc35b952780a19da5a38",
      "caae3f5bdac4ea0471e1e6905cf8ecf3fca847ff4fb9d0c9c9f4c961d09cbcff53797372"
      "f506a444424d0c61d3aa83f171e1d7907ec79e1a78c6b97a40fb6215"
      "ea26a4b18c6d28ce0e6062161a90cabad89777cf7a6ffbbcdfa8ceb33727c1f2e1c9c7f2"
      "b720fddc36e97c20b9466def419ecb3cb5f10e44a512b64bfa04f0e1",
      "0fd4f00e6f",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 8136U, &nameCaptureInfo);
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
  xInputs = emlrtCreateLogicalMatrix(1, 9);
  emlrtSetField(
      xEntryPoints, 0, "QualifiedName",
      emlrtMxCreateString("discretize_error_dynamics_FOH_RKV65_3DoF"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(9.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
          "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic "
          "SCP\\P"
          "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739774.07891203708));
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
                emlrtMxCreateString("nHjB0WTF6xfLPSgGXisb4C"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation
 * (_coder_discretize_error_dynamics_FOH_RKV65_3DoF_info.c) */
