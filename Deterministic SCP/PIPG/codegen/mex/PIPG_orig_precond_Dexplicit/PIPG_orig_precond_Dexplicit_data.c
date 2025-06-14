/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PIPG_orig_precond_Dexplicit_data.c
 *
 * Code generation for function 'PIPG_orig_precond_Dexplicit_data'
 *
 */

/* Include files */
#include "PIPG_orig_precond_Dexplicit_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;

const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;

emlrtContext emlrtContextGlobal = {
    true,                                                 /* bFirstTime */
    false,                                                /* bInitialized */
    131659U,                                              /* fVersionInfo */
    NULL,                                                 /* fErrorFunction */
    "PIPG_orig_precond_Dexplicit",                        /* fFunctionName */
    NULL,                                                 /* fRTCallStack */
    false,                                                /* bDebugMode */
    {2045744189U, 2170104910U, 2743257031U, 4284093946U}, /* fSigWrd */
    NULL                                                  /* fSigMem */
};

emlrtRSInfo k_emlrtRSI = {
    7,         /* lineNo */
    "getTime", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\+time\\getTime.m" /* pathName */
};

emlrtRSInfo l_emlrtRSI = {
    21,                     /* lineNo */
    "CoderTimeAPI/getTime", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\+time\\CoderTimeAPI.m" /* pathName */
};

emlrtRSInfo m_emlrtRSI = {
    148,                                  /* lineNo */
    "CoderTimeAPI/callEMLRTClockGettime", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\+time\\CoderTimeAPI.m" /* pathName */
};

emlrtRSInfo q_emlrtRSI = {
    35,              /* lineNo */
    "sparse/mtimes", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\mtimes.m" /* pathName */
};

emlrtRSInfo r_emlrtRSI = {
    216,        /* lineNo */
    "sfmtimes", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\mtimes.m" /* pathName */
};

emlrtRSInfo s_emlrtRSI = {
    241,               /* lineNo */
    "sfmtimes_vector", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\mtimes.m" /* pathName */
};

emlrtRSInfo t_emlrtRSI = {
    20,                               /* lineNo */
    "eml_int_forloop_overflow_check", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\eml\\eml_int_forloop_"
    "overflow_check.m" /* pathName */
};

emlrtRTEInfo c_emlrtRTEI = {
    159,                    /* lineNo */
    13,                     /* colNo */
    "coderTimeCheckStatus", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2024b\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\+time\\CoderTimeAPI.m" /* pName */
};

const char_T cv[26] = {'e', 'm', 'l', 'r', 't', 'C', 'l', 'o', 'c',
                       'k', 'G', 'e', 't', 't', 'i', 'm', 'e', 'M',
                       'o', 'n', 'o', 't', 'o', 'n', 'i', 'c'};

/* End of code generation (PIPG_orig_precond_Dexplicit_data.c) */
