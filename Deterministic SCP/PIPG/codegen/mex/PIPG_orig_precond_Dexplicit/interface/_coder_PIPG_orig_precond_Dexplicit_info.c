/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_PIPG_orig_precond_Dexplicit_info.c
 *
 * Code generation for function 'PIPG_orig_precond_Dexplicit'
 *
 */

/* Include files */
#include "_coder_PIPG_orig_precond_Dexplicit_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[7] = {
      "789ced984d8fd24018c707b3ba1b1395837af2c0176023fb62b2ded896854d5c6878590d"
      "5b03b51d9642db299d8105be80de4cbc78f0b277cf7e02133f8289f1"
      "e65efc087bb4d00e2f933425808565790e4c1ffed3f9cd3c347f9e14848e4f42008087c0"
      "89cb67cef8c0cdc3ee78078c07ab87dcf12e9383c1f71b63f751fda3",
      "3bcac820b04d9cc4907438b85341ba6a4806c9774c082c8891d6824a5fa9a81accab3acc"
      "8d26e95ea61f8d4883a427f5aeb92a94ebb9a60eac2a1eee501b4d06"
      "f5287b9c77c3a71e6cb0f560e7cdcaf3ab7fd8550e254debce81c7c6a4e7bb9e9247d7ff"
      "e4c3a3fa19f7967b291630b4b048aa12a9881903f296da82228fe4a6",
      "0e0d8245a169294d18e134096388c5783c11d93f781ee7e2bd4be72a9a2348ae4a98a872"
      "34c709d12c92eb90445f4986a21ae7220f09b4ece753ed4d88d81344"
      "e158488a1c325ab01db1077bd97ec9b7f5f13a943dcef968c23ab0e370fe567ffcfaf377"
      "2848de63fcfe2a481e8d45f1da1eeb4dfa1c3ff5e08519bd1dc3b1bd",
      "ccce453c7ff83adfde3dc5ddba1a4b0cf721f870fcf6013cf2a0d62f7bdcffbffcc79c92"
      "47d7bfefc3a37ace76070d126474c7f9e529f96c04e5b79f7d78545f"
      "22bf1d96de36dd55f5dbcd6f3bd741f268acbadf563276eff906c9d92394ba68d48ab546"
      "3696e4d67ecb4650fdeca627cf514c0bd5a04c16d63fd389d3f2be78",
      "f2c6f5b3ccf218acdbd08668ed01b83d3eb4eabca07cb6504c9fbe48e6f85a3ad528d4ce"
      "3b82b42bc8c9b5cfb2b1ac3ebbe83e7a56dfbdf4e1537d997c77b4b1"
      "a5e67b5b7c69d5796bdf9dcffa81f5638945da42efa3842cf5bc645a50468652e261dbd4"
      "545925ee8bc679fd2f6c31f9b01e8e8209324d7ba7f3e2ddf3e4398a",
      "829aef34383cdfdf19791f3c79e3fa627f6f5ae56d3d287ffaf523d8f71957dfff9c04c9"
      "a37153fdfe89072fcce83a91f9fd035c48c6f68ae96ca71ab3520a06"
      "37dfefff019982e471",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 7352U, &nameCaptureInfo);
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
  xInputs = emlrtCreateLogicalMatrix(1, 20);
  emlrtSetField(xEntryPoints, 0, "QualifiedName",
                emlrtMxCreateString("PIPG_orig_precond_Dexplicit"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(20.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString(
          "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
          "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic "
          "SCP\\P"
          "IPG\\PIPG_orig_precond_Dexplicit.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739774.10680555552));
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

/* End of code generation (_coder_PIPG_orig_precond_Dexplicit_info.c) */
