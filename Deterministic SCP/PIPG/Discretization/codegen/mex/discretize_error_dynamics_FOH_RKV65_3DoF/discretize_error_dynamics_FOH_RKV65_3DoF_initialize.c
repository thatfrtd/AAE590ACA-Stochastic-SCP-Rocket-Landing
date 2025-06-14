/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * discretize_error_dynamics_FOH_RKV65_3DoF_initialize.c
 *
 * Code generation for function
 * 'discretize_error_dynamics_FOH_RKV65_3DoF_initialize'
 *
 */

/* Include files */
#include "discretize_error_dynamics_FOH_RKV65_3DoF_initialize.h"
#include "_coder_discretize_error_dynamics_FOH_RKV65_3DoF_mex.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_data.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>

/* Function Declarations */
static void c_discretize_error_dynamics_FOH(void);

/* Function Definitions */
static void c_discretize_error_dynamics_FOH(void)
{
  static const int32_T lineInfo[84] = {
      9,   10,  11,  13,  14,  16,  17,  18,  19,  20,  21,  22,  23,  28,
      29,  31,  33,  35,  36,  38,  39,  47,  49,  51,  52,  53,  54,  55,
      57,  59,  60,  61,  62,  63,  64,  66,  68,  69,  71,  72,  74,  76,
      77,  78,  80,  81,  83,  84,  85,  86,  87,  89,  90,  91,  92,  93,
      94,  95,  97,  98,  100, 101, 102, 103, 104, 106, 107, 108, 109, 110,
      111, 112, 114, 115, 116, 117, 118, 119, 121, 122, 123, 124, 126, 127};
  static const int32_T b_lineInfo[14] = {130, 132, 133, 134, 135, 136, 137,
                                         138, 139, 141, 142, 143, 144, 145};
  static const int32_T b_offsets[4] = {0, 1, 2, 3};
  static const int32_T c_offsets[4] = {0, 1, 2, 3};
  static const int32_T d_offsets[4] = {0, 1, 2, 3};
  static const int32_T offsets[4] = {0, 1, 2, 3};
  __m128i r;
  __m128i r1;
  __m128i r2;
  __m128i r3;
  int32_T iv[19];
  int32_T iv2[14];
  int32_T iv1[12];
  int32_T iv4[4];
  int32_T iv3[2];
  mex_InitInfAndNan();
  zero_if_empty_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Helper Functions\\Ge"
      "neral\\zero_if_empty.m>zero_if_empty(codegen)";
  S_3DoF_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Dynamics Models\\3Do"
      "F\\S_3DoF.m>S_3DoF(codegen)";
  c_SymDynamics3DoF_mass_polar_co =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Dynamics Models\\3Do"
      "F\\SymDynamics3DoF_mass_polar.m>SymDynamics3DoF_mass_polar(codegen)";
  B_3DoF_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Dynamics Models\\3Do"
      "F\\B_3DoF.m>B_3DoF(codegen)";
  A_3DoF_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Dynamics Models\\3Do"
      "F\\A_3DoF.m>A_3DoF(codegen)";
  STM_diff_eq_FOH_complete_name =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m>STM_"
      "diff_eq_FOH(codegen)";
  d_discretize_error_dynamics_FOH =
      "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
      "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
      "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m>"
      "discretize_error_dynamics_FOH_RKV65_3DoF(codegen)";
  isMexOutdated = emlrtProfilerCheckMEXOutdated();
  emlrtProfilerRegisterMEXFcn(
      (char_T *)d_discretize_error_dynamics_FOH,
      (char_T
           *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
             "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic "
             "SCP\\P"
             "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",
      (char_T *)"discretize_error_dynamics_FOH_RKV65_3DoF", 84, &lineInfo[0],
      isMexOutdated);
  emlrtProfilerRegisterMEXFcn(
      (char_T *)STM_diff_eq_FOH_complete_name,
      (char_T
           *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
             "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic "
             "SCP\\P"
             "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",
      (char_T *)"discretize_error_dynamics_FOH_RKV65_3DoF>STM_diff_eq_FOH", 14,
      &b_lineInfo[0], isMexOutdated);
  r = _mm_set1_epi32(8);
  r1 = _mm_set1_epi32(0);
  r2 = _mm_loadu_si128((const __m128i *)&offsets[0]);
  _mm_storeu_si128((__m128i *)&iv[0], _mm_add_epi32(r, _mm_add_epi32(r1, r2)));
  r3 = _mm_set1_epi32(4);
  _mm_storeu_si128((__m128i *)&iv[4], _mm_add_epi32(r, _mm_add_epi32(r3, r2)));
  _mm_storeu_si128((__m128i *)&iv[8], _mm_add_epi32(r, _mm_add_epi32(r, r2)));
  _mm_storeu_si128((__m128i *)&iv[12],
                   _mm_add_epi32(r, _mm_add_epi32(_mm_set1_epi32(12), r2)));
  iv[16] = 24;
  iv[17] = 25;
  iv[18] = 26;
  emlrtProfilerRegisterMEXFcn(
      (char_T *)A_3DoF_complete_name,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Dynamics "
                "Models\\3Do"
                "F\\A_3DoF.m",
      (char_T *)"A_3DoF", 19, &iv[0], isMexOutdated);
  r2 = _mm_loadu_si128((const __m128i *)&b_offsets[0]);
  _mm_storeu_si128((__m128i *)&iv1[0], _mm_add_epi32(r, _mm_add_epi32(r1, r2)));
  _mm_storeu_si128((__m128i *)&iv1[4], _mm_add_epi32(r, _mm_add_epi32(r3, r2)));
  _mm_storeu_si128((__m128i *)&iv1[8], _mm_add_epi32(r, _mm_add_epi32(r, r2)));
  emlrtProfilerRegisterMEXFcn(
      (char_T *)B_3DoF_complete_name,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Dynamics "
                "Models\\3Do"
                "F\\B_3DoF.m",
      (char_T *)"B_3DoF", 12, &iv1[0], isMexOutdated);
  r2 = _mm_loadu_si128((const __m128i *)&c_offsets[0]);
  _mm_storeu_si128((__m128i *)&iv2[0], _mm_add_epi32(r, _mm_add_epi32(r1, r2)));
  _mm_storeu_si128((__m128i *)&iv2[4], _mm_add_epi32(r, _mm_add_epi32(r3, r2)));
  _mm_storeu_si128((__m128i *)&iv2[8], _mm_add_epi32(r, _mm_add_epi32(r, r2)));
  iv2[12] = 20;
  iv2[13] = 21;
  emlrtProfilerRegisterMEXFcn(
      (char_T *)c_SymDynamics3DoF_mass_polar_co,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Dynamics "
                "Models\\3Do"
                "F\\SymDynamics3DoF_mass_polar.m",
      (char_T *)"SymDynamics3DoF_mass_polar", 14, &iv2[0], isMexOutdated);
  iv3[0] = 8;
  iv3[1] = 9;
  emlrtProfilerRegisterMEXFcn(
      (char_T *)S_3DoF_complete_name,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Dynamics "
                "Models\\3Do"
                "F\\S_3DoF.m",
      (char_T *)"S_3DoF", 2, &iv3[0], isMexOutdated);
  _mm_storeu_si128(
      (__m128i *)&iv4[0],
      _mm_add_epi32(r3, _mm_add_epi32(r1, _mm_loadu_si128((
                                              const __m128i *)&d_offsets[0]))));
  emlrtProfilerRegisterMEXFcn(
      (char_T *)zero_if_empty_complete_name,
      (char_T *)"C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE "
                "590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Helper "
                "Functions\\Ge"
                "neral\\zero_if_empty.m",
      (char_T *)"zero_if_empty", 4, &iv4[0], isMexOutdated);
}

void discretize_error_dynamics_FOH_RKV65_3DoF_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    c_discretize_error_dynamics_FOH();
  }
  emlrtCheckProfilerStatus();
}

/* End of code generation
 * (discretize_error_dynamics_FOH_RKV65_3DoF_initialize.c) */
