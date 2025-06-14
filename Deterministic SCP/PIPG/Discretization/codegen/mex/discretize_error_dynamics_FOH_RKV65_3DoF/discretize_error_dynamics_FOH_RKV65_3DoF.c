/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * discretize_error_dynamics_FOH_RKV65_3DoF.c
 *
 * Code generation for function 'discretize_error_dynamics_FOH_RKV65_3DoF'
 *
 */

/* Include files */
#include "discretize_error_dynamics_FOH_RKV65_3DoF.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_data.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_emxutil.h"
#include "discretize_error_dynamics_FOH_RKV65_3DoF_types.h"
#include "eye.h"
#include "linspace.h"
#include "pagemtimes.h"
#include "reshapeSizeChecks.h"
#include "rt_nonfinite.h"
#include "vecnorm.h"
#include "zero_if_empty.h"
#include "mwmathutil.h"
#include <emmintrin.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 13,    /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 68,  /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 69,  /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 74,  /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 77,  /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 78,  /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 97,  /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 121, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 122, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 123, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 124, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 126, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 51,  /* lineNo */
  "reshapeSizeChecks",                 /* fcnName */
  "C:\\Program Files\\MATLAB\\R2024b\\toolbox\\eml\\eml\\+coder\\+internal\\reshapeSizeChecks.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 136, /* lineNo */
  "STM_diff_eq_FOH",                   /* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 141, /* lineNo */
  "STM_diff_eq_FOH",                   /* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 142, /* lineNo */
  "STM_diff_eq_FOH",                   /* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 143, /* lineNo */
  "STM_diff_eq_FOH",                   /* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 144, /* lineNo */
  "STM_diff_eq_FOH",                   /* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 80,  /* lineNo */
  13,                                  /* colNo */
  "reshapeSizeChecks",                 /* fName */
  "C:\\Program Files\\MATLAB\\R2024b\\toolbox\\eml\\eml\\+coder\\+internal\\reshapeSizeChecks.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 85,/* lineNo */
  23,                                  /* colNo */
  "reshapeSizeChecks",                 /* fName */
  "C:\\Program Files\\MATLAB\\R2024b\\toolbox\\eml\\eml\\+coder\\+internal\\reshapeSizeChecks.m"/* pName */
};

static emlrtRTEInfo c_emlrtRTEI = { 87,/* lineNo */
  23,                                  /* colNo */
  "reshapeSizeChecks",                 /* fName */
  "C:\\Program Files\\MATLAB\\R2024b\\toolbox\\eml\\eml\\+coder\\+internal\\reshapeSizeChecks.m"/* pName */
};

static emlrtECInfo emlrtECI = { 2,     /* nDims */
  126,                                 /* lineNo */
  21,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtDCInfo emlrtDCI = { 126,   /* lineNo */
  42,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtECInfo b_emlrtECI = { 3,   /* nDims */
  121,                                 /* lineNo */
  11,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo c_emlrtECI = { 1,   /* nDims */
  121,                                 /* lineNo */
  11,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo d_emlrtECI = { 3,   /* nDims */
  121,                                 /* lineNo */
  43,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo e_emlrtECI = { 2,   /* nDims */
  121,                                 /* lineNo */
  43,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo f_emlrtECI = { 1,   /* nDims */
  121,                                 /* lineNo */
  43,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtDCInfo b_emlrtDCI = { 123, /* lineNo */
  53,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo emlrtBCI = { 1,     /* iFirst */
  15,                                  /* iLast */
  122,                                 /* lineNo */
  55,                                  /* colNo */
  "u_ref",                             /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = { 122, /* lineNo */
  55,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { 1,   /* iFirst */
  15,                                  /* iLast */
  121,                                 /* lineNo */
  81,                                  /* colNo */
  "x_ref",                             /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = { 121, /* lineNo */
  81,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtECInfo g_emlrtECI = { 3,   /* nDims */
  111,                                 /* lineNo */
  21,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo h_emlrtECI = { 3,   /* nDims */
  110,                                 /* lineNo */
  27,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo i_emlrtECI = { 3,   /* nDims */
  109,                                 /* lineNo */
  26,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo j_emlrtECI = { 3,   /* nDims */
  108,                                 /* lineNo */
  21,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo k_emlrtECI = { 2,   /* nDims */
  107,                                 /* lineNo */
  21,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo l_emlrtECI = { -1,  /* nDims */
  97,                                  /* lineNo */
  112,                                 /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo m_emlrtECI = { -1,  /* nDims */
  97,                                  /* lineNo */
  83,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo n_emlrtECI = { -1,  /* nDims */
  97,                                  /* lineNo */
  55,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo o_emlrtECI = { -1,  /* nDims */
  97,                                  /* lineNo */
  32,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo p_emlrtECI = { -1,  /* nDims */
  97,                                  /* lineNo */
  14,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo q_emlrtECI = { 3,   /* nDims */
  94,                                  /* lineNo */
  27,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo r_emlrtECI = { -1,  /* nDims */
  74,                                  /* lineNo */
  108,                                 /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo s_emlrtECI = { 3,   /* nDims */
  93,                                  /* lineNo */
  33,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo t_emlrtECI = { -1,  /* nDims */
  74,                                  /* lineNo */
  79,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo u_emlrtECI = { 3,   /* nDims */
  92,                                  /* lineNo */
  32,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo v_emlrtECI = { -1,  /* nDims */
  74,                                  /* lineNo */
  51,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo w_emlrtECI = { 3,   /* nDims */
  91,                                  /* lineNo */
  27,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo x_emlrtECI = { -1,  /* nDims */
  74,                                  /* lineNo */
  28,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo y_emlrtECI = { 2,   /* nDims */
  90,                                  /* lineNo */
  25,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo ab_emlrtECI = { -1, /* nDims */
  74,                                  /* lineNo */
  10,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo bb_emlrtECI = { 2,  /* nDims */
  80,                                  /* lineNo */
  19,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo cb_emlrtECI = { 2,  /* nDims */
  81,                                  /* lineNo */
  19,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo db_emlrtECI = { 2,  /* nDims */
  81,                                  /* lineNo */
  20,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo eb_emlrtECI = { 2,  /* nDims */
  78,                                  /* lineNo */
  38,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo fb_emlrtECI = { 2,  /* nDims */
  71,                                  /* lineNo */
  15,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo gb_emlrtECI = { 2,  /* nDims */
  72,                                  /* lineNo */
  15,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo hb_emlrtECI = { 2,  /* nDims */
  72,                                  /* lineNo */
  16,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo ib_emlrtECI = { 2,  /* nDims */
  69,                                  /* lineNo */
  34,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 59,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtBCInfo c_emlrtBCI = { 1,   /* iFirst */
  15,                                  /* iLast */
  36,                                  /* lineNo */
  28,                                  /* colNo */
  "u_ref",                             /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = { 36,  /* lineNo */
  28,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { 1,   /* iFirst */
  15,                                  /* iLast */
  35,                                  /* lineNo */
  29,                                  /* colNo */
  "u_ref",                             /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo f_emlrtDCI = { 35,  /* lineNo */
  29,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { 1,   /* iFirst */
  15,                                  /* iLast */
  33,                                  /* lineNo */
  25,                                  /* colNo */
  "x_ref",                             /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo g_emlrtDCI = { 33,  /* lineNo */
  25,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  31,                                  /* lineNo */
  24,                                  /* colNo */
  "t",                                 /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  31,                                  /* lineNo */
  17,                                  /* colNo */
  "t",                                 /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  28,                                  /* lineNo */
  22,                                  /* colNo */
  "t",                                 /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  28,                                  /* lineNo */
  17,                                  /* colNo */
  "t",                                 /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo e_emlrtRTEI = { 17,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtBCInfo j_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  14,                                  /* lineNo */
  19,                                  /* colNo */
  "t",                                 /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo k_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  14,                                  /* lineNo */
  12,                                  /* colNo */
  "t",                                 /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo h_emlrtDCI = { 16,  /* lineNo */
  17,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo i_emlrtDCI = { 16,  /* lineNo */
  17,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo j_emlrtDCI = { 20,  /* lineNo */
  22,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  18,                                  /* lineNo */
  19,                                  /* colNo */
  "A_k",                               /* aName */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo k_emlrtDCI = { 21,  /* lineNo */
  23,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo l_emlrtDCI = { 22,  /* lineNo */
  17,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo m_emlrtDCI = { 23,  /* lineNo */
  17,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo n_emlrtDCI = { 51,  /* lineNo */
  21,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo o_emlrtDCI = { 52,  /* lineNo */
  23,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo p_emlrtDCI = { 53,  /* lineNo */
  28,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo q_emlrtDCI = { 54,  /* lineNo */
  29,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo r_emlrtDCI = { 55,  /* lineNo */
  23,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo s_emlrtDCI = { 90,  /* lineNo */
  27,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo t_emlrtDCI = { 91,  /* lineNo */
  29,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo u_emlrtDCI = { 92,  /* lineNo */
  34,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo v_emlrtDCI = { 93,  /* lineNo */
  35,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo w_emlrtDCI = { 94,  /* lineNo */
  29,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo m_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  136,                                 /* lineNo */
  33,                                  /* colNo */
  "t",                                 /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  136,                                 /* lineNo */
  42,                                  /* colNo */
  "x",                                 /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo o_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  136,                                 /* lineNo */
  51,                                  /* colNo */
  "u",                                 /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo p_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  137,                                 /* lineNo */
  33,                                  /* colNo */
  "t",                                 /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo q_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  137,                                 /* lineNo */
  42,                                  /* colNo */
  "x",                                 /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo r_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  137,                                 /* lineNo */
  51,                                  /* colNo */
  "u",                                 /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo s_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  138,                                 /* lineNo */
  51,                                  /* colNo */
  "t",                                 /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo t_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  138,                                 /* lineNo */
  60,                                  /* colNo */
  "x",                                 /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo u_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  138,                                 /* lineNo */
  69,                                  /* colNo */
  "u",                                 /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo jb_emlrtECI = { 3,  /* nDims */
  142,                                 /* lineNo */
  49,                                  /* colNo */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo kb_emlrtECI = { 1,  /* nDims */
  142,                                 /* lineNo */
  19,                                  /* colNo */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo lb_emlrtECI = { 2,  /* nDims */
  142,                                 /* lineNo */
  19,                                  /* colNo */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo mb_emlrtECI = { 3,  /* nDims */
  142,                                 /* lineNo */
  19,                                  /* colNo */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo nb_emlrtECI = { 3,  /* nDims */
  143,                                 /* lineNo */
  51,                                  /* colNo */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo ob_emlrtECI = { 1,  /* nDims */
  143,                                 /* lineNo */
  20,                                  /* colNo */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo pb_emlrtECI = { 2,  /* nDims */
  143,                                 /* lineNo */
  20,                                  /* colNo */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo qb_emlrtECI = { 3,  /* nDims */
  143,                                 /* lineNo */
  20,                                  /* colNo */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtECInfo rb_emlrtECI = { 1,  /* nDims */
  144,                                 /* lineNo */
  14,                                  /* colNo */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtBCInfo v_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  136,                                 /* lineNo */
  19,                                  /* colNo */
  "A_t",                               /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo w_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  137,                                 /* lineNo */
  19,                                  /* colNo */
  "B_t",                               /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo x_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  138,                                 /* lineNo */
  17,                                  /* colNo */
  "xdot",                              /* aName */
  "STM_diff_eq_FOH",                   /* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo l_emlrtRTEI = { 16,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo m_emlrtRTEI = { 20,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo n_emlrtRTEI = { 21,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo o_emlrtRTEI = { 22,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo p_emlrtRTEI = { 33,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo q_emlrtRTEI = { 51,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo r_emlrtRTEI = { 52,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo s_emlrtRTEI = { 53,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo t_emlrtRTEI = { 54,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo u_emlrtRTEI = { 55,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo v_emlrtRTEI = { 66,/* lineNo */
  9,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo w_emlrtRTEI = { 68,/* lineNo */
  32,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo x_emlrtRTEI = { 69,/* lineNo */
  34,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo y_emlrtRTEI = { 69,/* lineNo */
  33,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = { 121,/* lineNo */
  59,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = { 71,/* lineNo */
  15,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo cb_emlrtRTEI = { 72,/* lineNo */
  15,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo db_emlrtRTEI = { 74,/* lineNo */
  149,                                 /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = { 76,/* lineNo */
  13,                                  /* colNo */
  "eml_mtimes_helper",                 /* fName */
  "C:\\Program Files\\MATLAB\\R2024b\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_helper.m"/* pName */
};

static emlrtRTEInfo fb_emlrtRTEI = { 121,/* lineNo */
  43,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo gb_emlrtRTEI = { 77,/* lineNo */
  36,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo hb_emlrtRTEI = { 121,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ib_emlrtRTEI = { 107,/* lineNo */
  29,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo jb_emlrtRTEI = { 126,/* lineNo */
  21,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo kb_emlrtRTEI = { 107,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo lb_emlrtRTEI = { 108,/* lineNo */
  29,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo mb_emlrtRTEI = { 78,/* lineNo */
  37,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo nb_emlrtRTEI = { 108,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ob_emlrtRTEI = { 109,/* lineNo */
  39,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = { 109,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo qb_emlrtRTEI = { 110,/* lineNo */
  41,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo rb_emlrtRTEI = { 80,/* lineNo */
  19,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo sb_emlrtRTEI = { 81,/* lineNo */
  19,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo tb_emlrtRTEI = { 110,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ub_emlrtRTEI = { 83,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo vb_emlrtRTEI = { 111,/* lineNo */
  29,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo wb_emlrtRTEI = { 84,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo xb_emlrtRTEI = { 111,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo yb_emlrtRTEI = { 85,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ac_emlrtRTEI = { 86,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo bc_emlrtRTEI = { 87,/* lineNo */
  13,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo cc_emlrtRTEI = { 90,/* lineNo */
  25,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo dc_emlrtRTEI = { 97,/* lineNo */
  153,                                 /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ec_emlrtRTEI = { 90,/* lineNo */
  17,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo fc_emlrtRTEI = { 91,/* lineNo */
  27,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo gc_emlrtRTEI = { 91,/* lineNo */
  17,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo hc_emlrtRTEI = { 92,/* lineNo */
  32,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ic_emlrtRTEI = { 92,/* lineNo */
  17,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo jc_emlrtRTEI = { 93,/* lineNo */
  33,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo kc_emlrtRTEI = { 93,/* lineNo */
  17,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo lc_emlrtRTEI = { 94,/* lineNo */
  27,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo mc_emlrtRTEI = { 94,/* lineNo */
  17,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo nc_emlrtRTEI = { 13,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo oc_emlrtRTEI = { 60,/* lineNo */
  9,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo pc_emlrtRTEI = { 61,/* lineNo */
  9,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo qc_emlrtRTEI = { 62,/* lineNo */
  9,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo rc_emlrtRTEI = { 63,/* lineNo */
  9,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo sc_emlrtRTEI = { 64,/* lineNo */
  9,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo tc_emlrtRTEI = { 1,/* lineNo */
  56,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo vc_emlrtRTEI = { 132,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo wc_emlrtRTEI = { 133,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo xc_emlrtRTEI = { 134,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo yc_emlrtRTEI = { 142,/* lineNo */
  49,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ad_emlrtRTEI = { 142,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo bd_emlrtRTEI = { 143,/* lineNo */
  51,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo cd_emlrtRTEI = { 143,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo dd_emlrtRTEI = { 144,/* lineNo */
  5,                                   /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ed_emlrtRTEI = { 142,/* lineNo */
  19,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo jd_emlrtRTEI = { 111,/* lineNo */
  21,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo kd_emlrtRTEI = { 110,/* lineNo */
  27,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo ld_emlrtRTEI = { 108,/* lineNo */
  21,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo md_emlrtRTEI = { 107,/* lineNo */
  21,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRTEInfo nd_emlrtRTEI = { 143,/* lineNo */
  20,                                  /* colNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pName */
};

static emlrtRSInfo ib_emlrtRSI = { 94, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo jb_emlrtRSI = { 93, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo kb_emlrtRSI = { 92, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo lb_emlrtRSI = { 91, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo mb_emlrtRSI = { 90, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo nb_emlrtRSI = { 111,/* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo ob_emlrtRSI = { 110,/* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo pb_emlrtRSI = { 109,/* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo qb_emlrtRSI = { 108,/* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo rb_emlrtRSI = { 107,/* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo tb_emlrtRSI = { 80, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo ub_emlrtRSI = { 71, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo vb_emlrtRSI = { 81, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

static emlrtRSInfo wb_emlrtRSI = { 72, /* lineNo */
  "discretize_error_dynamics_FOH_RKV65_3DoF",/* fcnName */
  "C:\\Users\\thatf\\OneDrive\\Documents\\Purdue Classes\\AAE 590ACA\\AAE590ACA-Stochastic-SCP-Rocket-Landing\\Deterministic SCP\\P"
  "IPG\\Discretization\\discretize_error_dynamics_FOH_RKV65_3DoF.m"/* pathName */
};

/* Function Declarations */
static void STM_diff_eq_FOH(const emlrtStack *sp, const emxArray_real_T *t,
  const emxArray_real_T *x, const emxArray_real_T *u, const emxArray_real_T
  *sigma_plus, const emxArray_real_T *sigma_minus, const emxArray_real_T *STM,
  const emxArray_real_T *Phi_B_plus, const emxArray_real_T *Phi_B_minus, const
  emxArray_real_T *Phi_S, real_T L, real_T b_I, real_T alpha, emxArray_real_T
  *xdot, emxArray_real_T *A_kdot, emxArray_real_T *B_k_plusdot, emxArray_real_T *
  B_k_minusdot, emxArray_real_T *S_kdot);
static void b_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void b_times(const emlrtStack *sp, emxArray_real_T *in1, const
                    emxArray_real_T *in2, const emxArray_real_T *in3);
static void binary_expand_op(const emlrtStack *sp, const emlrtRSInfo in1, const
  emxArray_real_T *in2, const real_T in3[9], int32_T in4, real_T in5, const
  emxArray_real_T *in6, const emxArray_real_T *in7, const emxArray_real_T *in8,
  const emxArray_real_T *in9, const emxArray_real_T *in10, const emxArray_real_T
  *in11, const emxArray_real_T *in12, const emxArray_real_T *in13, const
  emxArray_real_T *in14, real_T in15, real_T in16, real_T in17, emxArray_real_T *
  in18, emxArray_real_T *in19, emxArray_real_T *in20, emxArray_real_T *in21,
  emxArray_real_T *in22);
static void binary_expand_op_1(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2, const real_T in3[30], int32_T in4, int32_T in5);
static void binary_expand_op_2(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2, const real_T in3[30], int32_T in4);
static void binary_expand_op_3(const emlrtStack *sp, const emlrtRSInfo in1,
  const emxArray_real_T *in2, real_T in3, const emxArray_real_T *in4, const
  emxArray_real_T *in5, const emxArray_real_T *in6, const emxArray_real_T *in7,
  const emxArray_real_T *in8, const emxArray_real_T *in9, const emxArray_real_T *
  in10, const emxArray_real_T *in11, const emxArray_real_T *in12, real_T in13,
  real_T in14, real_T in15, emxArray_real_T *in16, emxArray_real_T *in17,
  emxArray_real_T *in18, emxArray_real_T *in19, emxArray_real_T *in20);
static void binary_expand_op_6(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2, int32_T in3, int32_T in4, const emxArray_real_T *in5);
static void binary_expand_op_7(const emlrtStack *sp, emxArray_real_T *in1, const
  emlrtRSInfo in2, const emxArray_real_T *in3, const real_T in4[105], int32_T
  in5, int32_T in6);
static void binary_expand_op_8(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2, real_T in3, const emxArray_real_T *in4);
static void binary_expand_op_9(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2);
static void c_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void d_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void e_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void f_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void g_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void h_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void i_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void j_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void k_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2);
static void plus(const emlrtStack *sp, emxArray_real_T *in1, const
                 emxArray_real_T *in2);
static void times(const emlrtStack *sp, emxArray_real_T *in1, const
                  emxArray_real_T *in2);

/* Function Definitions */
static void STM_diff_eq_FOH(const emlrtStack *sp, const emxArray_real_T *t,
  const emxArray_real_T *x, const emxArray_real_T *u, const emxArray_real_T
  *sigma_plus, const emxArray_real_T *sigma_minus, const emxArray_real_T *STM,
  const emxArray_real_T *Phi_B_plus, const emxArray_real_T *Phi_B_minus, const
  emxArray_real_T *Phi_S, real_T L, real_T b_I, real_T alpha, emxArray_real_T
  *xdot, emxArray_real_T *A_kdot, emxArray_real_T *B_k_plusdot, emxArray_real_T *
  B_k_minusdot, emxArray_real_T *S_kdot)
{
  __m128d r1;
  __m128d r2;
  emlrtStack st;
  emxArray_real_T *A_t;
  emxArray_real_T *b_B_k_minusdot;
  emxArray_real_T *r;
  const real_T *sigma_minus_data;
  const real_T *sigma_plus_data;
  const real_T *u_data;
  const real_T *x_data;
  real_T *A_t_data;
  real_T *B_k_minusdot_data;
  real_T *xdot_data;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T k;
  int32_T loop_ub_tmp;
  int32_T t6_tmp_tmp;
  st.prev = sp;
  st.tls = sp->tls;
  sigma_minus_data = sigma_minus->data;
  sigma_plus_data = sigma_plus->data;
  u_data = u->data;
  x_data = x->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emlrtMEXProfilingFunctionEntry((char_T *)STM_diff_eq_FOH_complete_name,
    isMexOutdated);
  emlrtMEXProfilingStatement(1, isMexOutdated);
  emlrtMEXProfilingStatement(2, isMexOutdated);
  emxInit_real_T(sp, &A_t, 3, &vc_emlrtRTEI);
  i = A_t->size[0] * A_t->size[1] * A_t->size[2];
  A_t->size[0] = 7;
  A_t->size[1] = 7;
  A_t->size[2] = STM->size[2];
  emxEnsureCapacity_real_T(sp, A_t, i, &vc_emlrtRTEI);
  A_t_data = A_t->data;
  loop_ub_tmp = 49 * STM->size[2];
  for (i = 0; i < loop_ub_tmp; i++) {
    A_t_data[i] = 0.0;
  }

  emlrtMEXProfilingStatement(3, isMexOutdated);
  i = B_k_minusdot->size[0] * B_k_minusdot->size[1] * B_k_minusdot->size[2];
  B_k_minusdot->size[0] = 7;
  B_k_minusdot->size[1] = 2;
  B_k_minusdot->size[2] = Phi_B_plus->size[2];
  emxEnsureCapacity_real_T(sp, B_k_minusdot, i, &wc_emlrtRTEI);
  B_k_minusdot_data = B_k_minusdot->data;
  loop_ub_tmp = 14 * Phi_B_plus->size[2];
  for (i = 0; i < loop_ub_tmp; i++) {
    B_k_minusdot_data[i] = 0.0;
  }

  emlrtMEXProfilingStatement(4, isMexOutdated);
  i = xdot->size[0] * xdot->size[1];
  xdot->size[0] = 7;
  xdot->size[1] = x->size[1];
  emxEnsureCapacity_real_T(sp, xdot, i, &xc_emlrtRTEI);
  xdot_data = xdot->data;
  loop_ub_tmp = 7 * x->size[1];
  for (i = 0; i < loop_ub_tmp; i++) {
    xdot_data[i] = 0.0;
  }

  emlrtMEXProfilingStatement(5, isMexOutdated);
  i = t->size[1];
  for (k = 0; k < i; k++) {
    real_T A_t_tmp;
    real_T b_A_t_tmp;
    real_T t13_tmp;
    real_T t14;
    real_T t2_tmp;
    real_T t3_tmp;
    real_T t4_tmp;
    real_T t5_tmp;
    real_T t6;
    real_T t6_tmp;
    real_T t7;
    emlrtMEXProfilingStatement(6, isMexOutdated);
    if (k + 1 > i) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, i, &m_emlrtBCI, (emlrtConstCTX)sp);
    }

    st.site = &n_emlrtRSI;
    if (k + 1 > x->size[1]) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, x->size[1], &n_emlrtBCI, &st);
    }

    if (k + 1 > u->size[1]) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, u->size[1], &o_emlrtBCI, &st);
    }

    emlrtMEXProfilingFunctionEntry((char_T *)A_3DoF_complete_name, isMexOutdated);

    /* A_3DoF */
    /*     J_A = A_3DoF(T,IN2,IN3,L,I,ALPHA) */
    /*     This function was generated by the Symbolic Math Toolbox version 24.2. */
    /*     06-Jun-2025 01:33:53 */
    emlrtMEXProfilingStatement(1, isMexOutdated);
    emlrtMEXProfilingStatement(2, isMexOutdated);
    emlrtMEXProfilingStatement(3, isMexOutdated);
    emlrtMEXProfilingStatement(4, isMexOutdated);
    emlrtMEXProfilingStatement(5, isMexOutdated);
    loop_ub_tmp = 2 * k + 1;
    t2_tmp = muDoubleScalarCos(u_data[loop_ub_tmp]);
    emlrtMEXProfilingStatement(6, isMexOutdated);
    t3_tmp = muDoubleScalarSin(u_data[loop_ub_tmp]);
    emlrtMEXProfilingStatement(7, isMexOutdated);
    loop_ub_tmp = 7 * k + 4;
    t4_tmp = muDoubleScalarCos(x_data[loop_ub_tmp]);
    emlrtMEXProfilingStatement(8, isMexOutdated);
    t5_tmp = muDoubleScalarSin(x_data[loop_ub_tmp]);
    emlrtMEXProfilingStatement(9, isMexOutdated);
    t6_tmp_tmp = 7 * k + 6;
    t6_tmp = 1.0 / x_data[t6_tmp_tmp];
    emlrtMEXProfilingStatement(10, isMexOutdated);
    t7 = t6_tmp * t6_tmp;
    emlrtMEXProfilingStatement(11, isMexOutdated);
    emlrtMEXProfilingStatement(12, isMexOutdated);
    emlrtMEXProfilingStatement(13, isMexOutdated);
    emlrtMEXProfilingStatement(14, isMexOutdated);
    emlrtMEXProfilingStatement(15, isMexOutdated);
    emlrtMEXProfilingStatement(16, isMexOutdated);
    t6 = u_data[2 * k] * t2_tmp;
    t14 = u_data[2 * k] * t3_tmp;
    t13_tmp = t6 * t5_tmp + t14 * t4_tmp;
    emlrtMEXProfilingStatement(17, isMexOutdated);
    t14 = t6 * t4_tmp - t14 * t5_tmp;
    emlrtMEXProfilingStatement(18, isMexOutdated);
    if (k + 1 > A_t->size[2]) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, A_t->size[2], &v_emlrtBCI, &st);
    }

    A_t_data[49 * k] = 0.0;
    A_t_data[49 * k + 1] = 0.0;
    A_t_data[49 * k + 2] = 0.0;
    A_t_data[49 * k + 3] = 0.0;
    A_t_data[49 * k + 4] = 0.0;
    A_t_data[49 * k + 5] = 0.0;
    A_t_data[49 * k + 6] = 0.0;
    A_t_data[49 * k + 7] = 0.0;
    A_t_data[49 * k + 8] = 0.0;
    A_t_data[49 * k + 9] = 0.0;
    A_t_data[49 * k + 10] = 0.0;
    A_t_data[49 * k + 11] = 0.0;
    A_t_data[49 * k + 12] = 0.0;
    A_t_data[49 * k + 13] = 0.0;
    A_t_data[49 * k + 14] = 1.0;
    A_t_data[49 * k + 15] = 0.0;
    A_t_data[49 * k + 16] = 0.0;
    A_t_data[49 * k + 17] = 0.0;
    A_t_data[49 * k + 18] = 0.0;
    A_t_data[49 * k + 19] = 0.0;
    A_t_data[49 * k + 20] = 0.0;
    A_t_data[49 * k + 21] = 0.0;
    A_t_data[49 * k + 22] = 1.0;
    A_t_data[49 * k + 23] = 0.0;
    A_t_data[49 * k + 24] = 0.0;
    A_t_data[49 * k + 25] = 0.0;
    A_t_data[49 * k + 26] = 0.0;
    A_t_data[49 * k + 27] = 0.0;
    A_t_data[49 * k + 28] = 0.0;
    A_t_data[49 * k + 29] = 0.0;
    A_t_tmp = -t6_tmp * t13_tmp;
    A_t_data[49 * k + 30] = A_t_tmp;
    b_A_t_tmp = t6_tmp * t14;
    A_t_data[49 * k + 31] = b_A_t_tmp;
    A_t_data[49 * k + 32] = 0.0;
    A_t_data[49 * k + 33] = 0.0;
    A_t_data[49 * k + 34] = 0.0;
    A_t_data[49 * k + 35] = 0.0;
    A_t_data[49 * k + 36] = 0.0;
    A_t_data[49 * k + 37] = 0.0;
    A_t_data[49 * k + 38] = 0.0;
    A_t_data[49 * k + 39] = 1.0;
    A_t_data[49 * k + 40] = 0.0;
    A_t_data[49 * k + 41] = 0.0;
    A_t_data[49 * k + 42] = 0.0;
    A_t_data[49 * k + 43] = 0.0;
    A_t_data[49 * k + 44] = -t7 * t14;
    A_t_data[49 * k + 45] = -t7 * t13_tmp;
    A_t_data[49 * k + 46] = 0.0;
    A_t_data[49 * k + 47] = 0.0;
    A_t_data[49 * k + 48] = 0.0;
    emlrtMEXProfilingStatement(19, isMexOutdated);
    emlrtMEXProfilingFunctionExit(isMexOutdated);
    emlrtMEXProfilingStatement(7, isMexOutdated);
    if (k + 1 > i) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, i, &p_emlrtBCI, (emlrtConstCTX)sp);
    }

    if (k + 1 > x->size[1]) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, x->size[1], &q_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    if (k + 1 > u->size[1]) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, u->size[1], &r_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    emlrtMEXProfilingFunctionEntry((char_T *)B_3DoF_complete_name, isMexOutdated);

    /* B_3DoF */
    /*     J_B = B_3DoF(T,IN2,IN3,L,I,ALPHA) */
    /*     This function was generated by the Symbolic Math Toolbox version 24.2. */
    /*     06-Jun-2025 01:33:53 */
    emlrtMEXProfilingStatement(1, isMexOutdated);
    emlrtMEXProfilingStatement(2, isMexOutdated);
    emlrtMEXProfilingStatement(3, isMexOutdated);
    emlrtMEXProfilingStatement(4, isMexOutdated);
    emlrtMEXProfilingStatement(5, isMexOutdated);
    emlrtMEXProfilingStatement(6, isMexOutdated);
    emlrtMEXProfilingStatement(7, isMexOutdated);
    emlrtMEXProfilingStatement(8, isMexOutdated);
    emlrtMEXProfilingStatement(9, isMexOutdated);
    t6 = 1.0 / b_I;
    emlrtMEXProfilingStatement(10, isMexOutdated);
    emlrtMEXProfilingStatement(11, isMexOutdated);
    if (k + 1 > B_k_minusdot->size[2]) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, B_k_minusdot->size[2], &w_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    B_k_minusdot_data[14 * k] = 0.0;
    B_k_minusdot_data[14 * k + 1] = 0.0;
    B_k_minusdot_data[14 * k + 2] = t6_tmp * (t2_tmp * t4_tmp - t3_tmp * t5_tmp);
    B_k_minusdot_data[14 * k + 3] = t6_tmp * (t2_tmp * t5_tmp + t3_tmp * t4_tmp);
    B_k_minusdot_data[14 * k + 4] = 0.0;
    B_k_minusdot_data[14 * k + 5] = -L * t3_tmp * t6;
    B_k_minusdot_data[14 * k + 6] = -alpha;
    B_k_minusdot_data[14 * k + 7] = 0.0;
    B_k_minusdot_data[14 * k + 8] = 0.0;
    B_k_minusdot_data[14 * k + 9] = A_t_tmp;
    B_k_minusdot_data[14 * k + 10] = b_A_t_tmp;
    B_k_minusdot_data[14 * k + 11] = 0.0;
    B_k_minusdot_data[14 * k + 12] = -L * u_data[2 * k] * t2_tmp * t6;
    B_k_minusdot_data[14 * k + 13] = 0.0;
    emlrtMEXProfilingStatement(12, isMexOutdated);
    emlrtMEXProfilingFunctionExit(isMexOutdated);
    emlrtMEXProfilingStatement(8, isMexOutdated);
    if (k + 1 > i) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, i, &s_emlrtBCI, (emlrtConstCTX)sp);
    }

    if (k + 1 > x->size[1]) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, x->size[1], &t_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    if (k + 1 > u->size[1]) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, u->size[1], &u_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    emlrtMEXProfilingFunctionEntry((char_T *)c_SymDynamics3DoF_mass_polar_co,
      isMexOutdated);

    /* SymDynamics3DoF_mass_polar */
    /*     XDOT = SymDynamics3DoF_mass_polar(T,IN2,IN3,L,I,ALPHA) */
    /*     This function was generated by the Symbolic Math Toolbox version 24.2. */
    /*     06-Jun-2025 01:33:53 */
    emlrtMEXProfilingStatement(1, isMexOutdated);
    emlrtMEXProfilingStatement(2, isMexOutdated);
    emlrtMEXProfilingStatement(3, isMexOutdated);
    emlrtMEXProfilingStatement(4, isMexOutdated);
    emlrtMEXProfilingStatement(5, isMexOutdated);
    emlrtMEXProfilingStatement(6, isMexOutdated);
    emlrtMEXProfilingStatement(7, isMexOutdated);
    emlrtMEXProfilingStatement(8, isMexOutdated);
    emlrtMEXProfilingStatement(9, isMexOutdated);
    emlrtMEXProfilingStatement(10, isMexOutdated);
    emlrtMEXProfilingStatement(11, isMexOutdated);
    emlrtMEXProfilingStatement(12, isMexOutdated);
    emlrtMEXProfilingStatement(13, isMexOutdated);
    if (k + 1 > xdot->size[1]) {
      emlrtDynamicBoundsCheckR2012b(k + 1, 1, xdot->size[1], &x_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    i1 = 7 * k + 2;
    xdot_data[7 * k] = x_data[i1];
    i2 = 7 * k + 3;
    xdot_data[7 * k + 1] = x_data[i2];
    xdot_data[i1] = b_A_t_tmp;
    xdot_data[i2] = t6_tmp * t13_tmp - 0.0037114;
    i1 = 7 * k + 5;
    xdot_data[loop_ub_tmp] = x_data[i1];
    xdot_data[i1] = -(L * u_data[2 * k] * t3_tmp) / b_I;
    xdot_data[t6_tmp_tmp] = -u_data[2 * k] * alpha;
    emlrtMEXProfilingStatement(14, isMexOutdated);
    emlrtMEXProfilingFunctionExit(isMexOutdated);
    emlrtMEXProfilingStatement(9, isMexOutdated);
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }

  emlrtMEXProfilingStatement(10, isMexOutdated);
  st.site = &o_emlrtRSI;
  pagemtimes(&st, A_t, STM, A_kdot);
  emlrtMEXProfilingStatement(11, isMexOutdated);
  k = B_k_minusdot->size[2];
  if ((B_k_minusdot->size[2] != sigma_plus->size[2]) && ((B_k_minusdot->size[2]
        != 1) && (sigma_plus->size[2] != 1))) {
    emlrtDimSizeImpxCheckR2021b(B_k_minusdot->size[2], sigma_plus->size[2],
      &jb_emlrtECI, (emlrtConstCTX)sp);
  }

  emxInit_real_T(sp, &r, 3, &ed_emlrtRTEI);
  st.site = &p_emlrtRSI;
  b_pagemtimes(&st, A_t, Phi_B_plus, r);
  xdot_data = r->data;
  if ((r->size[0] != 7) && (r->size[0] != 1)) {
    emlrtDimSizeImpxCheckR2021b(r->size[0], 7, &kb_emlrtECI, (emlrtConstCTX)sp);
  }

  if ((r->size[1] != 2) && (r->size[1] != 1)) {
    emlrtDimSizeImpxCheckR2021b(r->size[1], 2, &lb_emlrtECI, (emlrtConstCTX)sp);
  }

  if (B_k_minusdot->size[2] == sigma_plus->size[2]) {
    i = B_k_plusdot->size[0] * B_k_plusdot->size[1] * B_k_plusdot->size[2];
    B_k_plusdot->size[0] = 7;
    B_k_plusdot->size[1] = 2;
    B_k_plusdot->size[2] = B_k_minusdot->size[2];
    emxEnsureCapacity_real_T(sp, B_k_plusdot, i, &yc_emlrtRTEI);
    A_t_data = B_k_plusdot->data;
    for (i = 0; i < k; i++) {
      for (i1 = 0; i1 < 2; i1++) {
        i2 = 7 * i1 + 14 * i;
        r1 = _mm_loadu_pd(&B_k_minusdot_data[i2]);
        r2 = _mm_set1_pd(sigma_plus_data[i]);
        _mm_storeu_pd(&A_t_data[i2], _mm_mul_pd(r1, r2));
        r1 = _mm_loadu_pd(&B_k_minusdot_data[i2 + 2]);
        _mm_storeu_pd(&A_t_data[i2 + 2], _mm_mul_pd(r1, r2));
        r1 = _mm_loadu_pd(&B_k_minusdot_data[i2 + 4]);
        _mm_storeu_pd(&A_t_data[i2 + 4], _mm_mul_pd(r1, r2));
        A_t_data[i2 + 6] = sigma_plus_data[i] * B_k_minusdot_data[i2 + 6];
      }
    }
  } else {
    st.site = &p_emlrtRSI;
    b_times(&st, B_k_plusdot, B_k_minusdot, sigma_plus);
  }

  if ((r->size[2] != B_k_plusdot->size[2]) && ((r->size[2] != 1) &&
       (B_k_plusdot->size[2] != 1))) {
    emlrtDimSizeImpxCheckR2021b(r->size[2], B_k_plusdot->size[2], &mb_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if (r->size[2] == B_k_plusdot->size[2]) {
    k = r->size[0] * r->size[1] * r->size[2];
    i = B_k_plusdot->size[0] * B_k_plusdot->size[1] * B_k_plusdot->size[2];
    B_k_plusdot->size[0] = r->size[0];
    B_k_plusdot->size[1] = r->size[1];
    B_k_plusdot->size[2] = r->size[2];
    emxEnsureCapacity_real_T(sp, B_k_plusdot, i, &ad_emlrtRTEI);
    A_t_data = B_k_plusdot->data;
    loop_ub_tmp = (k / 2) << 1;
    t6_tmp_tmp = loop_ub_tmp - 2;
    for (i = 0; i <= t6_tmp_tmp; i += 2) {
      r1 = _mm_loadu_pd(&xdot_data[i]);
      r2 = _mm_loadu_pd(&A_t_data[i]);
      _mm_storeu_pd(&A_t_data[i], _mm_add_pd(r1, r2));
    }

    for (i = loop_ub_tmp; i < k; i++) {
      A_t_data[i] += xdot_data[i];
    }
  } else {
    st.site = &p_emlrtRSI;
    k_plus(&st, B_k_plusdot, r);
  }

  emlrtMEXProfilingStatement(12, isMexOutdated);
  k = B_k_minusdot->size[2];
  if ((B_k_minusdot->size[2] != sigma_minus->size[2]) && ((B_k_minusdot->size[2]
        != 1) && (sigma_minus->size[2] != 1))) {
    emlrtDimSizeImpxCheckR2021b(B_k_minusdot->size[2], sigma_minus->size[2],
      &nb_emlrtECI, (emlrtConstCTX)sp);
  }

  st.site = &q_emlrtRSI;
  b_pagemtimes(&st, A_t, Phi_B_minus, r);
  xdot_data = r->data;
  if ((r->size[0] != 7) && (r->size[0] != 1)) {
    emlrtDimSizeImpxCheckR2021b(r->size[0], 7, &ob_emlrtECI, (emlrtConstCTX)sp);
  }

  if ((r->size[1] != 2) && (r->size[1] != 1)) {
    emlrtDimSizeImpxCheckR2021b(r->size[1], 2, &pb_emlrtECI, (emlrtConstCTX)sp);
  }

  if (B_k_minusdot->size[2] == sigma_minus->size[2]) {
    emxInit_real_T(sp, &b_B_k_minusdot, 3, &bd_emlrtRTEI);
    i = b_B_k_minusdot->size[0] * b_B_k_minusdot->size[1] * b_B_k_minusdot->
      size[2];
    b_B_k_minusdot->size[0] = 7;
    b_B_k_minusdot->size[1] = 2;
    b_B_k_minusdot->size[2] = B_k_minusdot->size[2];
    emxEnsureCapacity_real_T(sp, b_B_k_minusdot, i, &bd_emlrtRTEI);
    A_t_data = b_B_k_minusdot->data;
    for (i = 0; i < k; i++) {
      for (i1 = 0; i1 < 2; i1++) {
        i2 = 7 * i1 + 14 * i;
        r1 = _mm_loadu_pd(&B_k_minusdot_data[i2]);
        r2 = _mm_set1_pd(sigma_minus_data[i]);
        _mm_storeu_pd(&A_t_data[i2], _mm_mul_pd(r1, r2));
        r1 = _mm_loadu_pd(&B_k_minusdot_data[i2 + 2]);
        _mm_storeu_pd(&A_t_data[i2 + 2], _mm_mul_pd(r1, r2));
        r1 = _mm_loadu_pd(&B_k_minusdot_data[i2 + 4]);
        _mm_storeu_pd(&A_t_data[i2 + 4], _mm_mul_pd(r1, r2));
        A_t_data[i2 + 6] = sigma_minus_data[i] * B_k_minusdot_data[i2 + 6];
      }
    }

    i = B_k_minusdot->size[0] * B_k_minusdot->size[1] * B_k_minusdot->size[2];
    B_k_minusdot->size[0] = 7;
    B_k_minusdot->size[1] = 2;
    B_k_minusdot->size[2] = b_B_k_minusdot->size[2];
    emxEnsureCapacity_real_T(sp, B_k_minusdot, i, &bd_emlrtRTEI);
    B_k_minusdot_data = B_k_minusdot->data;
    loop_ub_tmp = 14 * b_B_k_minusdot->size[2];
    for (i = 0; i < loop_ub_tmp; i++) {
      B_k_minusdot_data[i] = A_t_data[i];
    }

    emxFree_real_T(sp, &b_B_k_minusdot);
  } else {
    st.site = &q_emlrtRSI;
    times(&st, B_k_minusdot, sigma_minus);
  }

  if ((r->size[2] != B_k_minusdot->size[2]) && ((r->size[2] != 1) &&
       (B_k_minusdot->size[2] != 1))) {
    emlrtDimSizeImpxCheckR2021b(r->size[2], B_k_minusdot->size[2], &qb_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if (r->size[2] == B_k_minusdot->size[2]) {
    k = r->size[0] * r->size[1] * r->size[2];
    i = B_k_minusdot->size[0] * B_k_minusdot->size[1] * B_k_minusdot->size[2];
    B_k_minusdot->size[0] = r->size[0];
    B_k_minusdot->size[1] = r->size[1];
    B_k_minusdot->size[2] = r->size[2];
    emxEnsureCapacity_real_T(sp, B_k_minusdot, i, &cd_emlrtRTEI);
    B_k_minusdot_data = B_k_minusdot->data;
    loop_ub_tmp = (k / 2) << 1;
    t6_tmp_tmp = loop_ub_tmp - 2;
    for (i = 0; i <= t6_tmp_tmp; i += 2) {
      r1 = _mm_loadu_pd(&xdot_data[i]);
      r2 = _mm_loadu_pd(&B_k_minusdot_data[i]);
      _mm_storeu_pd(&B_k_minusdot_data[i], _mm_add_pd(r1, r2));
    }

    for (i = loop_ub_tmp; i < k; i++) {
      B_k_minusdot_data[i] += xdot_data[i];
    }
  } else {
    st.site = &q_emlrtRSI;
    k_plus(&st, B_k_minusdot, r);
  }

  emlrtMEXProfilingStatement(13, isMexOutdated);
  st.site = &r_emlrtRSI;
  c_pagemtimes(&st, A_t, Phi_S, r);
  xdot_data = r->data;
  emxFree_real_T(sp, &A_t);
  if ((r->size[0] != 7) && (r->size[0] != 1)) {
    emlrtDimSizeImpxCheckR2021b(r->size[0], 7, &rb_emlrtECI, (emlrtConstCTX)sp);
  }

  if (r->size[0] == 7) {
    i = S_kdot->size[0] * S_kdot->size[1] * S_kdot->size[2];
    S_kdot->size[0] = 7;
    k = r->size[1];
    S_kdot->size[1] = r->size[1];
    loop_ub_tmp = r->size[2];
    S_kdot->size[2] = r->size[2];
    emxEnsureCapacity_real_T(sp, S_kdot, i, &dd_emlrtRTEI);
    A_t_data = S_kdot->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      for (i1 = 0; i1 < k; i1++) {
        for (i2 = 0; i2 < 7; i2++) {
          A_t_data[(i2 + 7 * i1) + 7 * S_kdot->size[1] * i] = xdot_data[(i2 +
            r->size[0] * i1) + r->size[0] * r->size[1] * i];
        }
      }
    }
  } else {
    st.site = &r_emlrtRSI;
    binary_expand_op_9(&st, S_kdot, r);
  }

  emxFree_real_T(sp, &r);
  emlrtMEXProfilingStatement(14, isMexOutdated);
  emlrtMEXProfilingFunctionExit(isMexOutdated);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void b_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in2;
  const real_T *in2_data;
  real_T *b_in2_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_2;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in2, 3, &jc_emlrtRTEI);
  i = b_in2->size[0] * b_in2->size[1] * b_in2->size[2];
  b_in2->size[0] = 7;
  b_in2->size[1] = 2;
  if (in1->size[2] == 1) {
    loop_ub = in2->size[2];
  } else {
    loop_ub = in1->size[2];
  }

  b_in2->size[2] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &jc_emlrtRTEI);
  b_in2_data = b_in2->data;
  stride_0_2 = (in2->size[2] != 1);
  stride_1_2 = (in1->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 2; i1++) {
      __m128d r;
      __m128d r1;
      int32_T i4;
      i2 = 7 * i1 + 14 * aux_0_2;
      r = _mm_loadu_pd(&in2_data[i2]);
      i3 = 7 * i1 + 14 * aux_1_2;
      r1 = _mm_loadu_pd(&in1_data[i3]);
      i4 = 7 * i1 + 14 * i;
      _mm_storeu_pd(&b_in2_data[i4], _mm_add_pd(r, r1));
      r = _mm_loadu_pd(&in2_data[i2 + 2]);
      r1 = _mm_loadu_pd(&in1_data[i3 + 2]);
      _mm_storeu_pd(&b_in2_data[i4 + 2], _mm_add_pd(r, r1));
      r = _mm_loadu_pd(&in2_data[i2 + 4]);
      r1 = _mm_loadu_pd(&in1_data[i3 + 4]);
      _mm_storeu_pd(&b_in2_data[i4 + 4], _mm_add_pd(r, r1));
      b_in2_data[i4 + 6] = in2_data[i2 + 6] + in1_data[i3 + 6];
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  in1->size[1] = 2;
  emxEnsureCapacity_real_T(sp, in1, i, &jc_emlrtRTEI);
  loop_ub = b_in2->size[2];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_in2->size[2];
  emxEnsureCapacity_real_T(sp, in1, i, &jc_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 2; i1++) {
      for (i2 = 0; i2 < 7; i2++) {
        i3 = (i2 + 7 * i1) + 14 * i;
        in1_data[i3] = b_in2_data[i3];
      }
    }
  }

  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void b_times(const emlrtStack *sp, emxArray_real_T *in1, const
                    emxArray_real_T *in2, const emxArray_real_T *in3)
{
  const real_T *in2_data;
  const real_T *in3_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T i;
  int32_T i1;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_2;
  in3_data = in3->data;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  in1->size[1] = 2;
  emxEnsureCapacity_real_T(sp, in1, i, &yc_emlrtRTEI);
  if (in3->size[2] == 1) {
    loop_ub = in2->size[2];
  } else {
    loop_ub = in3->size[2];
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &yc_emlrtRTEI);
  in1_data = in1->data;
  stride_0_2 = (in2->size[2] != 1);
  stride_1_2 = (in3->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 2; i1++) {
      __m128d r;
      int32_T i2;
      int32_T i3;
      i2 = 7 * i1 + 14 * aux_0_2;
      r = _mm_loadu_pd(&in2_data[i2]);
      i3 = 7 * i1 + 14 * i;
      _mm_storeu_pd(&in1_data[i3], _mm_mul_pd(r, _mm_set1_pd(in3_data[aux_1_2])));
      r = _mm_loadu_pd(&in2_data[i2 + 2]);
      _mm_storeu_pd(&in1_data[i3 + 2], _mm_mul_pd(r, _mm_set1_pd
        (in3_data[aux_1_2])));
      r = _mm_loadu_pd(&in2_data[i2 + 4]);
      _mm_storeu_pd(&in1_data[i3 + 4], _mm_mul_pd(r, _mm_set1_pd
        (in3_data[aux_1_2])));
      in1_data[i3 + 6] = in3_data[aux_1_2] * in2_data[i2 + 6];
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }
}

static void binary_expand_op(const emlrtStack *sp, const emlrtRSInfo in1, const
  emxArray_real_T *in2, const real_T in3[9], int32_T in4, real_T in5, const
  emxArray_real_T *in6, const emxArray_real_T *in7, const emxArray_real_T *in8,
  const emxArray_real_T *in9, const emxArray_real_T *in10, const emxArray_real_T
  *in11, const emxArray_real_T *in12, const emxArray_real_T *in13, const
  emxArray_real_T *in14, real_T in15, real_T in16, real_T in17, emxArray_real_T *
  in18, emxArray_real_T *in19, emxArray_real_T *in20, emxArray_real_T *in21,
  emxArray_real_T *in22)
{
  __m128d r;
  emlrtStack st;
  emxArray_real_T b_in10;
  emxArray_real_T b_in9;
  emxArray_real_T *b_in2;
  emxArray_real_T *b_in7;
  const real_T *in2_data;
  const real_T *in7_data;
  const real_T *in8_data;
  real_T b_in3;
  real_T *b_in2_data;
  int32_T iv[3];
  int32_T iv1[3];
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T b_unnamed_idx_2;
  int32_T i;
  int32_T loop_ub;
  int32_T scalarLB;
  int32_T unnamed_idx_2;
  int32_T vectorUB;
  st.prev = sp;
  st.tls = sp->tls;
  in8_data = in8->data;
  in7_data = in7->data;
  in2_data = in2->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  b_in3 = in3[in4 + 1] * in5;
  unnamed_idx_2 = in9->size[1];
  b_unnamed_idx_2 = in10->size[1];
  emxInit_real_T(sp, &b_in2, 2, &dc_emlrtRTEI);
  i = b_in2->size[0] * b_in2->size[1];
  b_in2->size[0] = 1;
  loop_ub = in2->size[1];
  b_in2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &dc_emlrtRTEI);
  b_in2_data = b_in2->data;
  scalarLB = (loop_ub / 2) << 1;
  vectorUB = scalarLB - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&in2_data[i]);
    _mm_storeu_pd(&b_in2_data[i], _mm_add_pd(r, _mm_set1_pd(b_in3)));
  }

  for (i = scalarLB; i < loop_ub; i++) {
    b_in2_data[i] = in2_data[i] + b_in3;
  }

  emxInit_real_T(sp, &b_in7, 2, &rb_emlrtRTEI);
  i = b_in7->size[0] * b_in7->size[1];
  b_in7->size[0] = 2;
  if (in8->size[1] == 1) {
    loop_ub = in7->size[1];
  } else {
    loop_ub = in8->size[1];
  }

  b_in7->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in7, i, &rb_emlrtRTEI);
  b_in2_data = b_in7->data;
  scalarLB = (in7->size[1] != 1);
  vectorUB = (in8->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < loop_ub; i++) {
    __m128d r1;
    r = _mm_loadu_pd(&in7_data[2 * aux_0_1]);
    r1 = _mm_loadu_pd(&in8_data[2 * aux_1_1]);
    _mm_storeu_pd(&b_in2_data[2 * i], _mm_add_pd(r, r1));
    aux_1_1 += vectorUB;
    aux_0_1 += scalarLB;
  }

  b_in9 = *in9;
  iv[0] = 1;
  iv[1] = 1;
  iv[2] = unnamed_idx_2;
  b_in9.size = &iv[0];
  b_in9.numDimensions = 3;
  b_in10 = *in10;
  iv1[0] = 1;
  iv1[1] = 1;
  iv1[2] = b_unnamed_idx_2;
  b_in10.size = &iv1[0];
  b_in10.numDimensions = 3;
  st.site = (emlrtRSInfo *)&in1;
  STM_diff_eq_FOH(&st, b_in2, in6, b_in7, &b_in9, &b_in10, in11, in12, in13,
                  in14, in15, in16, in17, in18, in19, in20, in21, in22);
  emxFree_real_T(sp, &b_in7);
  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void binary_expand_op_1(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2, const real_T in3[30], int32_T in4, int32_T in5)
{
  const real_T *in2_data;
  real_T *in1_data;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 2;
  emxEnsureCapacity_real_T(sp, in1, i, &sb_emlrtRTEI);
  i = in5 - in4;
  if (i == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = i;
  }

  stride_0_1 = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, stride_0_1, &sb_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (i != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < loop_ub; i++) {
    __m128d r;
    r = _mm_loadu_pd(&in3[(in4 + aux_1_1) << 1]);
    _mm_storeu_pd(&in1_data[2 * i], _mm_mul_pd(_mm_set1_pd(in2_data[aux_0_1]), r));
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

static void binary_expand_op_2(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2, const real_T in3[30], int32_T in4)
{
  const real_T *in2_data;
  real_T *in1_data;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 2;
  emxEnsureCapacity_real_T(sp, in1, i, &rb_emlrtRTEI);
  if (in4 == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in4;
  }

  i = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &rb_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in4 != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < loop_ub; i++) {
    __m128d r;
    r = _mm_loadu_pd(&in3[aux_1_1 << 1]);
    _mm_storeu_pd(&in1_data[2 * i], _mm_mul_pd(_mm_set1_pd(in2_data[aux_0_1]), r));
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

static void binary_expand_op_3(const emlrtStack *sp, const emlrtRSInfo in1,
  const emxArray_real_T *in2, real_T in3, const emxArray_real_T *in4, const
  emxArray_real_T *in5, const emxArray_real_T *in6, const emxArray_real_T *in7,
  const emxArray_real_T *in8, const emxArray_real_T *in9, const emxArray_real_T *
  in10, const emxArray_real_T *in11, const emxArray_real_T *in12, real_T in13,
  real_T in14, real_T in15, emxArray_real_T *in16, emxArray_real_T *in17,
  emxArray_real_T *in18, emxArray_real_T *in19, emxArray_real_T *in20)
{
  __m128d r;
  emlrtStack st;
  emxArray_real_T b_in7;
  emxArray_real_T b_in8;
  emxArray_real_T *b_in2;
  emxArray_real_T *b_in5;
  const real_T *in2_data;
  const real_T *in5_data;
  const real_T *in6_data;
  real_T d;
  real_T *b_in2_data;
  int32_T iv[3];
  int32_T iv1[3];
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T b_unnamed_idx_2;
  int32_T i;
  int32_T loop_ub;
  int32_T scalarLB;
  int32_T unnamed_idx_2;
  int32_T vectorUB;
  st.prev = sp;
  st.tls = sp->tls;
  in6_data = in6->data;
  in5_data = in5->data;
  in2_data = in2->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  d = 0.0 * in3;
  unnamed_idx_2 = in7->size[1];
  b_unnamed_idx_2 = in8->size[1];
  emxInit_real_T(sp, &b_in2, 2, &db_emlrtRTEI);
  i = b_in2->size[0] * b_in2->size[1];
  b_in2->size[0] = 1;
  loop_ub = in2->size[1];
  b_in2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &db_emlrtRTEI);
  b_in2_data = b_in2->data;
  scalarLB = (loop_ub / 2) << 1;
  vectorUB = scalarLB - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&in2_data[i]);
    _mm_storeu_pd(&b_in2_data[i], _mm_add_pd(r, _mm_set1_pd(d)));
  }

  for (i = scalarLB; i < loop_ub; i++) {
    b_in2_data[i] = in2_data[i] + d;
  }

  emxInit_real_T(sp, &b_in5, 2, &bb_emlrtRTEI);
  i = b_in5->size[0] * b_in5->size[1];
  b_in5->size[0] = 2;
  if (in6->size[1] == 1) {
    loop_ub = in5->size[1];
  } else {
    loop_ub = in6->size[1];
  }

  b_in5->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in5, i, &bb_emlrtRTEI);
  b_in2_data = b_in5->data;
  scalarLB = (in5->size[1] != 1);
  vectorUB = (in6->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < loop_ub; i++) {
    __m128d r1;
    r = _mm_loadu_pd(&in5_data[2 * aux_0_1]);
    r1 = _mm_loadu_pd(&in6_data[2 * aux_1_1]);
    _mm_storeu_pd(&b_in2_data[2 * i], _mm_add_pd(r, r1));
    aux_1_1 += vectorUB;
    aux_0_1 += scalarLB;
  }

  b_in7 = *in7;
  iv[0] = 1;
  iv[1] = 1;
  iv[2] = unnamed_idx_2;
  b_in7.size = &iv[0];
  b_in7.numDimensions = 3;
  b_in8 = *in8;
  iv1[0] = 1;
  iv1[1] = 1;
  iv1[2] = b_unnamed_idx_2;
  b_in8.size = &iv1[0];
  b_in8.numDimensions = 3;
  st.site = (emlrtRSInfo *)&in1;
  STM_diff_eq_FOH(&st, b_in2, in4, b_in5, &b_in7, &b_in8, in9, in10, in11, in12,
                  in13, in14, in15, in16, in17, in18, in19, in20);
  emxFree_real_T(sp, &b_in5);
  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void binary_expand_op_6(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2, int32_T in3, int32_T in4, const emxArray_real_T *in5)
{
  const real_T *in2_data;
  const real_T *in5_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in5_data = in5->data;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  emxEnsureCapacity_real_T(sp, in1, i, &x_emlrtRTEI);
  i = in4 - in3;
  if (in5->size[1] == 1) {
    loop_ub = i;
  } else {
    loop_ub = in5->size[1];
  }

  stride_0_1 = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, stride_0_1, &x_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (i != 1);
  stride_1_1 = (in5->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = in2_data[in3 + i * stride_0_1] - in5_data[i * stride_1_1];
  }
}

static void binary_expand_op_7(const emlrtStack *sp, emxArray_real_T *in1, const
  emlrtRSInfo in2, const emxArray_real_T *in3, const real_T in4[105], int32_T
  in5, int32_T in6)
{
  emlrtStack st;
  emxArray_real_T *b_in3;
  const real_T *in3_data;
  real_T *b_in3_data;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T i;
  int32_T i1;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  st.prev = sp;
  st.tls = sp->tls;
  in3_data = in3->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in3, 2, &jb_emlrtRTEI);
  i = b_in3->size[0] * b_in3->size[1];
  b_in3->size[0] = 7;
  i1 = in6 - in5;
  if (i1 == 1) {
    loop_ub = in3->size[1];
  } else {
    loop_ub = i1;
  }

  b_in3->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in3, i, &jb_emlrtRTEI);
  b_in3_data = b_in3->data;
  stride_0_1 = (in3->size[1] != 1);
  stride_1_1 = (i1 != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < loop_ub; i++) {
    __m128d r;
    __m128d r1;
    r = _mm_loadu_pd(&in3_data[7 * aux_0_1]);
    i1 = 7 * (in5 + aux_1_1);
    r1 = _mm_loadu_pd(&in4[i1]);
    _mm_storeu_pd(&b_in3_data[7 * i], _mm_sub_pd(r, r1));
    r = _mm_loadu_pd(&in3_data[7 * aux_0_1 + 2]);
    r1 = _mm_loadu_pd(&in4[i1 + 2]);
    _mm_storeu_pd(&b_in3_data[7 * i + 2], _mm_sub_pd(r, r1));
    r = _mm_loadu_pd(&in3_data[7 * aux_0_1 + 4]);
    r1 = _mm_loadu_pd(&in4[i1 + 4]);
    _mm_storeu_pd(&b_in3_data[7 * i + 4], _mm_sub_pd(r, r1));
    b_in3_data[7 * i + 6] = in3_data[7 * aux_0_1 + 6] - in4[i1 + 6];
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }

  st.site = (emlrtRSInfo *)&in2;
  vecnorm(&st, b_in3, in1);
  emxFree_real_T(sp, &b_in3);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void binary_expand_op_8(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2, real_T in3, const emxArray_real_T *in4)
{
  const real_T *in2_data;
  const real_T *in4_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T b_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_0;
  int32_T stride_1_2;
  in4_data = in4->data;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  emxEnsureCapacity_real_T(sp, in1, i, &hb_emlrtRTEI);
  loop_ub = in4->size[1];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &hb_emlrtRTEI);
  if (in4->size[2] == 1) {
    b_loop_ub = (int32_T)(in3 - 1.0);
  } else {
    b_loop_ub = in4->size[2];
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &hb_emlrtRTEI);
  in1_data = in1->data;
  stride_0_2 = ((int32_T)(in3 - 1.0) != 1);
  stride_1_0 = (in4->size[0] != 1);
  stride_1_2 = (in4->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      for (i2 = 0; i2 < 7; i2++) {
        in1_data[(i2 + 7 * i1) + 7 * in1->size[1] * i] = in2_data[i2 + 7 *
          aux_0_2] - in4_data[(i2 * stride_1_0 + in4->size[0] * i1) + in4->size
          [0] * in4->size[1] * aux_1_2];
      }
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }
}

static void binary_expand_op_9(const emlrtStack *sp, emxArray_real_T *in1, const
  emxArray_real_T *in2)
{
  const real_T *in2_data;
  real_T *in1_data;
  int32_T b_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T loop_ub;
  int32_T stride_0_0;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  emxEnsureCapacity_real_T(sp, in1, i, &dd_emlrtRTEI);
  loop_ub = in2->size[1];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &dd_emlrtRTEI);
  b_loop_ub = in2->size[2];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &dd_emlrtRTEI);
  in1_data = in1->data;
  stride_0_0 = (in2->size[0] != 1);
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      for (i2 = 0; i2 < 7; i2++) {
        in1_data[(i2 + 7 * i1) + 7 * in1->size[1] * i] = in2_data[(i2 *
          stride_0_0 + in2->size[0] * i1) + in2->size[0] * in2->size[1] * i];
      }
    }
  }
}

static void c_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in2;
  const real_T *in2_data;
  real_T *b_in2_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_2;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in2, 3, &fc_emlrtRTEI);
  i = b_in2->size[0] * b_in2->size[1] * b_in2->size[2];
  b_in2->size[0] = 7;
  b_in2->size[1] = 7;
  if (in1->size[2] == 1) {
    loop_ub = in2->size[2];
  } else {
    loop_ub = in1->size[2];
  }

  b_in2->size[2] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &fc_emlrtRTEI);
  b_in2_data = b_in2->data;
  stride_0_2 = (in2->size[2] != 1);
  stride_1_2 = (in1->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 7; i1++) {
      __m128d r;
      __m128d r1;
      int32_T i4;
      i2 = 7 * i1 + 49 * aux_0_2;
      r = _mm_loadu_pd(&in2_data[i2]);
      i3 = 7 * i1 + 49 * aux_1_2;
      r1 = _mm_loadu_pd(&in1_data[i3]);
      i4 = 7 * i1 + 49 * i;
      _mm_storeu_pd(&b_in2_data[i4], _mm_add_pd(r, r1));
      r = _mm_loadu_pd(&in2_data[i2 + 2]);
      r1 = _mm_loadu_pd(&in1_data[i3 + 2]);
      _mm_storeu_pd(&b_in2_data[i4 + 2], _mm_add_pd(r, r1));
      r = _mm_loadu_pd(&in2_data[i2 + 4]);
      r1 = _mm_loadu_pd(&in1_data[i3 + 4]);
      _mm_storeu_pd(&b_in2_data[i4 + 4], _mm_add_pd(r, r1));
      b_in2_data[i4 + 6] = in2_data[i2 + 6] + in1_data[i3 + 6];
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  in1->size[1] = 7;
  emxEnsureCapacity_real_T(sp, in1, i, &fc_emlrtRTEI);
  loop_ub = b_in2->size[2];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_in2->size[2];
  emxEnsureCapacity_real_T(sp, in1, i, &fc_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 7; i1++) {
      for (i2 = 0; i2 < 7; i2++) {
        i3 = (i2 + 7 * i1) + 49 * i;
        in1_data[i3] = b_in2_data[i3];
      }
    }
  }

  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void d_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in2;
  const real_T *in2_data;
  real_T *b_in2_data;
  real_T *in1_data;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in2, 2, &cc_emlrtRTEI);
  i = b_in2->size[0] * b_in2->size[1];
  b_in2->size[0] = 7;
  if (in1->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in1->size[1];
  }

  b_in2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &cc_emlrtRTEI);
  b_in2_data = b_in2->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in1->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < loop_ub; i++) {
    __m128d r;
    __m128d r1;
    r = _mm_loadu_pd(&in2_data[7 * aux_0_1]);
    r1 = _mm_loadu_pd(&in1_data[7 * aux_1_1]);
    _mm_storeu_pd(&b_in2_data[7 * i], _mm_add_pd(r, r1));
    r = _mm_loadu_pd(&in2_data[7 * aux_0_1 + 2]);
    r1 = _mm_loadu_pd(&in1_data[7 * aux_1_1 + 2]);
    _mm_storeu_pd(&b_in2_data[7 * i + 2], _mm_add_pd(r, r1));
    r = _mm_loadu_pd(&in2_data[7 * aux_0_1 + 4]);
    r1 = _mm_loadu_pd(&in1_data[7 * aux_1_1 + 4]);
    _mm_storeu_pd(&b_in2_data[7 * i + 4], _mm_add_pd(r, r1));
    b_in2_data[7 * i + 6] = in2_data[7 * aux_0_1 + 6] + in1_data[7 * aux_1_1 + 6];
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }

  i = in1->size[0] * in1->size[1];
  in1->size[0] = 7;
  emxEnsureCapacity_real_T(sp, in1, i, &cc_emlrtRTEI);
  loop_ub = b_in2->size[1];
  i = in1->size[0] * in1->size[1];
  in1->size[1] = b_in2->size[1];
  emxEnsureCapacity_real_T(sp, in1, i, &cc_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (stride_1_1 = 0; stride_1_1 < 7; stride_1_1++) {
      stride_0_1 = stride_1_1 + 7 * i;
      in1_data[stride_0_1] = b_in2_data[stride_0_1];
    }
  }

  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void e_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_2;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, 3, &jd_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1] * b_in1->size[2];
  b_in1->size[0] = 7;
  b_in1->size[1] = 1;
  if (in2->size[2] == 1) {
    loop_ub = in1->size[2];
  } else {
    loop_ub = in2->size[2];
  }

  b_in1->size[2] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &jd_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_2 = (in1->size[2] != 1);
  stride_1_2 = (in2->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < loop_ub; i++) {
    __m128d r;
    __m128d r1;
    r = _mm_loadu_pd(&in1_data[7 * aux_0_2]);
    r1 = _mm_loadu_pd(&in2_data[7 * aux_1_2]);
    _mm_storeu_pd(&b_in1_data[7 * i], _mm_add_pd(r, r1));
    r = _mm_loadu_pd(&in1_data[7 * aux_0_2 + 2]);
    r1 = _mm_loadu_pd(&in2_data[7 * aux_1_2 + 2]);
    _mm_storeu_pd(&b_in1_data[7 * i + 2], _mm_add_pd(r, r1));
    r = _mm_loadu_pd(&in1_data[7 * aux_0_2 + 4]);
    r1 = _mm_loadu_pd(&in2_data[7 * aux_1_2 + 4]);
    _mm_storeu_pd(&b_in1_data[7 * i + 4], _mm_add_pd(r, r1));
    b_in1_data[7 * i + 6] = in1_data[7 * aux_0_2 + 6] + in2_data[7 * aux_1_2 + 6];
    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  in1->size[1] = 1;
  emxEnsureCapacity_real_T(sp, in1, i, &jd_emlrtRTEI);
  loop_ub = b_in1->size[2];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_in1->size[2];
  emxEnsureCapacity_real_T(sp, in1, i, &jd_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (stride_1_2 = 0; stride_1_2 < 7; stride_1_2++) {
      stride_0_2 = stride_1_2 + 7 * i;
      in1_data[stride_0_2] = b_in1_data[stride_0_2];
    }
  }

  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void f_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_2;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, 3, &kd_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1] * b_in1->size[2];
  b_in1->size[0] = 7;
  b_in1->size[1] = 2;
  if (in2->size[2] == 1) {
    loop_ub = in1->size[2];
  } else {
    loop_ub = in2->size[2];
  }

  b_in1->size[2] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &kd_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_2 = (in1->size[2] != 1);
  stride_1_2 = (in2->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 2; i1++) {
      __m128d r;
      __m128d r1;
      int32_T i4;
      i2 = 7 * i1 + 14 * aux_0_2;
      r = _mm_loadu_pd(&in1_data[i2]);
      i3 = 7 * i1 + 14 * aux_1_2;
      r1 = _mm_loadu_pd(&in2_data[i3]);
      i4 = 7 * i1 + 14 * i;
      _mm_storeu_pd(&b_in1_data[i4], _mm_add_pd(r, r1));
      r = _mm_loadu_pd(&in1_data[i2 + 2]);
      r1 = _mm_loadu_pd(&in2_data[i3 + 2]);
      _mm_storeu_pd(&b_in1_data[i4 + 2], _mm_add_pd(r, r1));
      r = _mm_loadu_pd(&in1_data[i2 + 4]);
      r1 = _mm_loadu_pd(&in2_data[i3 + 4]);
      _mm_storeu_pd(&b_in1_data[i4 + 4], _mm_add_pd(r, r1));
      b_in1_data[i4 + 6] = in1_data[i2 + 6] + in2_data[i3 + 6];
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  in1->size[1] = 2;
  emxEnsureCapacity_real_T(sp, in1, i, &kd_emlrtRTEI);
  loop_ub = b_in1->size[2];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_in1->size[2];
  emxEnsureCapacity_real_T(sp, in1, i, &kd_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 2; i1++) {
      for (i2 = 0; i2 < 7; i2++) {
        i3 = (i2 + 7 * i1) + 14 * i;
        in1_data[i3] = b_in1_data[i3];
      }
    }
  }

  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void g_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_2;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, 3, &ld_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1] * b_in1->size[2];
  b_in1->size[0] = 7;
  b_in1->size[1] = 7;
  if (in2->size[2] == 1) {
    loop_ub = in1->size[2];
  } else {
    loop_ub = in2->size[2];
  }

  b_in1->size[2] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &ld_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_2 = (in1->size[2] != 1);
  stride_1_2 = (in2->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 7; i1++) {
      __m128d r;
      __m128d r1;
      int32_T i4;
      i2 = 7 * i1 + 49 * aux_0_2;
      r = _mm_loadu_pd(&in1_data[i2]);
      i3 = 7 * i1 + 49 * aux_1_2;
      r1 = _mm_loadu_pd(&in2_data[i3]);
      i4 = 7 * i1 + 49 * i;
      _mm_storeu_pd(&b_in1_data[i4], _mm_add_pd(r, r1));
      r = _mm_loadu_pd(&in1_data[i2 + 2]);
      r1 = _mm_loadu_pd(&in2_data[i3 + 2]);
      _mm_storeu_pd(&b_in1_data[i4 + 2], _mm_add_pd(r, r1));
      r = _mm_loadu_pd(&in1_data[i2 + 4]);
      r1 = _mm_loadu_pd(&in2_data[i3 + 4]);
      _mm_storeu_pd(&b_in1_data[i4 + 4], _mm_add_pd(r, r1));
      b_in1_data[i4 + 6] = in1_data[i2 + 6] + in2_data[i3 + 6];
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  in1->size[1] = 7;
  emxEnsureCapacity_real_T(sp, in1, i, &ld_emlrtRTEI);
  loop_ub = b_in1->size[2];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_in1->size[2];
  emxEnsureCapacity_real_T(sp, in1, i, &ld_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 7; i1++) {
      for (i2 = 0; i2 < 7; i2++) {
        i3 = (i2 + 7 * i1) + 49 * i;
        in1_data[i3] = b_in1_data[i3];
      }
    }
  }

  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void h_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, 2, &md_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 7;
  if (in2->size[1] == 1) {
    loop_ub = in1->size[1];
  } else {
    loop_ub = in2->size[1];
  }

  b_in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &md_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in2->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < loop_ub; i++) {
    __m128d r;
    __m128d r1;
    r = _mm_loadu_pd(&in1_data[7 * aux_0_1]);
    r1 = _mm_loadu_pd(&in2_data[7 * aux_1_1]);
    _mm_storeu_pd(&b_in1_data[7 * i], _mm_add_pd(r, r1));
    r = _mm_loadu_pd(&in1_data[7 * aux_0_1 + 2]);
    r1 = _mm_loadu_pd(&in2_data[7 * aux_1_1 + 2]);
    _mm_storeu_pd(&b_in1_data[7 * i + 2], _mm_add_pd(r, r1));
    r = _mm_loadu_pd(&in1_data[7 * aux_0_1 + 4]);
    r1 = _mm_loadu_pd(&in2_data[7 * aux_1_1 + 4]);
    _mm_storeu_pd(&b_in1_data[7 * i + 4], _mm_add_pd(r, r1));
    b_in1_data[7 * i + 6] = in1_data[7 * aux_0_1 + 6] + in2_data[7 * aux_1_1 + 6];
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }

  i = in1->size[0] * in1->size[1];
  in1->size[0] = 7;
  emxEnsureCapacity_real_T(sp, in1, i, &md_emlrtRTEI);
  loop_ub = b_in1->size[1];
  i = in1->size[0] * in1->size[1];
  in1->size[1] = b_in1->size[1];
  emxEnsureCapacity_real_T(sp, in1, i, &md_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (stride_1_1 = 0; stride_1_1 < 7; stride_1_1++) {
      stride_0_1 = stride_1_1 + 7 * i;
      in1_data[stride_0_1] = b_in1_data[stride_0_1];
    }
  }

  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void i_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in2;
  const real_T *in2_data;
  real_T *b_in2_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T b_loop_ub;
  int32_T c_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T loop_ub;
  int32_T stride_0_0;
  int32_T stride_0_1;
  int32_T stride_0_2;
  int32_T stride_1_0;
  int32_T stride_1_1;
  int32_T stride_1_2;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in2, 3, &fb_emlrtRTEI);
  if (in1->size[0] == 1) {
    loop_ub = in2->size[0];
  } else {
    loop_ub = in1->size[0];
  }

  i = b_in2->size[0] * b_in2->size[1] * b_in2->size[2];
  b_in2->size[0] = loop_ub;
  if (in1->size[1] == 1) {
    b_loop_ub = in2->size[1];
  } else {
    b_loop_ub = in1->size[1];
  }

  b_in2->size[1] = b_loop_ub;
  if (in1->size[2] == 1) {
    c_loop_ub = in2->size[2];
  } else {
    c_loop_ub = in1->size[2];
  }

  b_in2->size[2] = c_loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &fb_emlrtRTEI);
  b_in2_data = b_in2->data;
  stride_0_0 = (in2->size[0] != 1);
  stride_0_1 = (in2->size[1] != 1);
  stride_0_2 = (in2->size[2] != 1);
  stride_1_0 = (in1->size[0] != 1);
  stride_1_1 = (in1->size[1] != 1);
  stride_1_2 = (in1->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < c_loop_ub; i++) {
    int32_T aux_0_1;
    int32_T aux_1_1;
    aux_0_1 = 0;
    aux_1_1 = 0;
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_in2_data[(i2 + b_in2->size[0] * i1) + b_in2->size[0] * b_in2->size[1] *
          i] = in2_data[(i2 * stride_0_0 + in2->size[0] * aux_0_1) + in2->size[0]
          * in2->size[1] * aux_0_2] + in1_data[(i2 * stride_1_0 + in1->size[0] *
          aux_1_1) + in1->size[0] * in1->size[1] * aux_1_2];
      }

      aux_1_1 += stride_1_1;
      aux_0_1 += stride_0_1;
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = loop_ub;
  in1->size[1] = b_loop_ub;
  in1->size[2] = c_loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &fb_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < c_loop_ub; i++) {
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        in1_data[(i2 + in1->size[0] * i1) + in1->size[0] * in1->size[1] * i] =
          b_in2_data[(i2 + b_in2->size[0] * i1) + b_in2->size[0] * b_in2->size[1]
          * i];
      }
    }
  }

  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void j_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T b_loop_ub;
  int32_T c_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T loop_ub;
  int32_T stride_0_0;
  int32_T stride_0_1;
  int32_T stride_0_2;
  int32_T stride_1_0;
  int32_T stride_1_1;
  int32_T stride_1_2;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, 3, &fb_emlrtRTEI);
  if (in2->size[0] == 1) {
    loop_ub = in1->size[0];
  } else {
    loop_ub = in2->size[0];
  }

  i = b_in1->size[0] * b_in1->size[1] * b_in1->size[2];
  b_in1->size[0] = loop_ub;
  if (in2->size[1] == 1) {
    b_loop_ub = in1->size[1];
  } else {
    b_loop_ub = in2->size[1];
  }

  b_in1->size[1] = b_loop_ub;
  if (in2->size[2] == 1) {
    c_loop_ub = in1->size[2];
  } else {
    c_loop_ub = in2->size[2];
  }

  b_in1->size[2] = c_loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &fb_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_0_1 = (in1->size[1] != 1);
  stride_0_2 = (in1->size[2] != 1);
  stride_1_0 = (in2->size[0] != 1);
  stride_1_1 = (in2->size[1] != 1);
  stride_1_2 = (in2->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < c_loop_ub; i++) {
    int32_T aux_0_1;
    int32_T aux_1_1;
    aux_0_1 = 0;
    aux_1_1 = 0;
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_in1_data[(i2 + b_in1->size[0] * i1) + b_in1->size[0] * b_in1->size[1] *
          i] = in1_data[(i2 * stride_0_0 + in1->size[0] * aux_0_1) + in1->size[0]
          * in1->size[1] * aux_0_2] + in2_data[(i2 * stride_1_0 + in2->size[0] *
          aux_1_1) + in2->size[0] * in2->size[1] * aux_1_2];
      }

      aux_1_1 += stride_1_1;
      aux_0_1 += stride_0_1;
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = loop_ub;
  in1->size[1] = b_loop_ub;
  in1->size[2] = c_loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &fb_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < c_loop_ub; i++) {
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        in1_data[(i2 + in1->size[0] * i1) + in1->size[0] * in1->size[1] * i] =
          b_in1_data[(i2 + b_in1->size[0] * i1) + b_in1->size[0] * b_in1->size[1]
          * i];
      }
    }
  }

  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void k_plus(const emlrtStack *sp, emxArray_real_T *in1, const
                   emxArray_real_T *in2)
{
  emxArray_real_T *b_in2;
  const real_T *in2_data;
  real_T *b_in2_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T b_loop_ub;
  int32_T c_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_2;
  int32_T vectorUB;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in2, 3, &nd_emlrtRTEI);
  loop_ub = in2->size[0];
  i = b_in2->size[0] * b_in2->size[1] * b_in2->size[2];
  b_in2->size[0] = loop_ub;
  b_loop_ub = in2->size[1];
  b_in2->size[1] = b_loop_ub;
  if (in1->size[2] == 1) {
    c_loop_ub = in2->size[2];
  } else {
    c_loop_ub = in1->size[2];
  }

  b_in2->size[2] = c_loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &nd_emlrtRTEI);
  b_in2_data = b_in2->data;
  stride_0_2 = (in2->size[2] != 1);
  stride_1_2 = (in1->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < c_loop_ub; i++) {
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      int32_T scalarLB;
      scalarLB = (loop_ub / 2) << 1;
      vectorUB = scalarLB - 2;
      for (i2 = 0; i2 <= vectorUB; i2 += 2) {
        __m128d r;
        __m128d r1;
        r = _mm_loadu_pd(&in2_data[(i2 + in2->size[0] * i1) + in2->size[0] *
                         in2->size[1] * aux_0_2]);
        r1 = _mm_loadu_pd(&in1_data[(i2 + 7 * i1) + 14 * aux_1_2]);
        _mm_storeu_pd(&b_in2_data[(i2 + 7 * i1) + 14 * i], _mm_add_pd(r, r1));
      }

      for (i2 = scalarLB; i2 < loop_ub; i2++) {
        vectorUB = i2 + 7 * i1;
        b_in2_data[vectorUB + 14 * i] = in2_data[(i2 + in2->size[0] * i1) +
          in2->size[0] * in2->size[1] * aux_0_2] + in1_data[vectorUB + 14 *
          aux_1_2];
      }
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  in1->size[1] = 2;
  emxEnsureCapacity_real_T(sp, in1, i, &nd_emlrtRTEI);
  loop_ub = b_in2->size[2];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_in2->size[2];
  emxEnsureCapacity_real_T(sp, in1, i, &nd_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 2; i1++) {
      for (i2 = 0; i2 < 7; i2++) {
        vectorUB = (i2 + 7 * i1) + 14 * i;
        in1_data[vectorUB] = b_in2_data[vectorUB];
      }
    }
  }

  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void plus(const emlrtStack *sp, emxArray_real_T *in1, const
                 emxArray_real_T *in2)
{
  emxArray_real_T *b_in2;
  const real_T *in2_data;
  real_T *b_in2_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_2;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in2, 3, &lc_emlrtRTEI);
  i = b_in2->size[0] * b_in2->size[1] * b_in2->size[2];
  b_in2->size[0] = 7;
  b_in2->size[1] = 1;
  if (in1->size[2] == 1) {
    loop_ub = in2->size[2];
  } else {
    loop_ub = in1->size[2];
  }

  b_in2->size[2] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &lc_emlrtRTEI);
  b_in2_data = b_in2->data;
  stride_0_2 = (in2->size[2] != 1);
  stride_1_2 = (in1->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < loop_ub; i++) {
    __m128d r;
    __m128d r1;
    r = _mm_loadu_pd(&in2_data[7 * aux_0_2]);
    r1 = _mm_loadu_pd(&in1_data[7 * aux_1_2]);
    _mm_storeu_pd(&b_in2_data[7 * i], _mm_add_pd(r, r1));
    r = _mm_loadu_pd(&in2_data[7 * aux_0_2 + 2]);
    r1 = _mm_loadu_pd(&in1_data[7 * aux_1_2 + 2]);
    _mm_storeu_pd(&b_in2_data[7 * i + 2], _mm_add_pd(r, r1));
    r = _mm_loadu_pd(&in2_data[7 * aux_0_2 + 4]);
    r1 = _mm_loadu_pd(&in1_data[7 * aux_1_2 + 4]);
    _mm_storeu_pd(&b_in2_data[7 * i + 4], _mm_add_pd(r, r1));
    b_in2_data[7 * i + 6] = in2_data[7 * aux_0_2 + 6] + in1_data[7 * aux_1_2 + 6];
    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  in1->size[1] = 1;
  emxEnsureCapacity_real_T(sp, in1, i, &lc_emlrtRTEI);
  loop_ub = b_in2->size[2];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_in2->size[2];
  emxEnsureCapacity_real_T(sp, in1, i, &lc_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (stride_1_2 = 0; stride_1_2 < 7; stride_1_2++) {
      stride_0_2 = stride_1_2 + 7 * i;
      in1_data[stride_0_2] = b_in2_data[stride_0_2];
    }
  }

  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void times(const emlrtStack *sp, emxArray_real_T *in1, const
                  emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T aux_0_2;
  int32_T aux_1_2;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T loop_ub;
  int32_T stride_0_2;
  int32_T stride_1_2;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, 3, &bd_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1] * b_in1->size[2];
  b_in1->size[0] = 7;
  b_in1->size[1] = 2;
  if (in2->size[2] == 1) {
    loop_ub = in1->size[2];
  } else {
    loop_ub = in2->size[2];
  }

  b_in1->size[2] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &bd_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_2 = (in1->size[2] != 1);
  stride_1_2 = (in2->size[2] != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 2; i1++) {
      __m128d r;
      i2 = 7 * i1 + 14 * aux_0_2;
      r = _mm_loadu_pd(&in1_data[i2]);
      i3 = 7 * i1 + 14 * i;
      _mm_storeu_pd(&b_in1_data[i3], _mm_mul_pd(r, _mm_set1_pd(in2_data[aux_1_2])));
      r = _mm_loadu_pd(&in1_data[i2 + 2]);
      _mm_storeu_pd(&b_in1_data[i3 + 2], _mm_mul_pd(r, _mm_set1_pd
        (in2_data[aux_1_2])));
      r = _mm_loadu_pd(&in1_data[i2 + 4]);
      _mm_storeu_pd(&b_in1_data[i3 + 4], _mm_mul_pd(r, _mm_set1_pd
        (in2_data[aux_1_2])));
      b_in1_data[i3 + 6] = in2_data[aux_1_2] * in1_data[i2 + 6];
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[0] = 7;
  in1->size[1] = 2;
  emxEnsureCapacity_real_T(sp, in1, i, &bd_emlrtRTEI);
  loop_ub = b_in1->size[2];
  i = in1->size[0] * in1->size[1] * in1->size[2];
  in1->size[2] = b_in1->size[2];
  emxEnsureCapacity_real_T(sp, in1, i, &bd_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 2; i1++) {
      for (i2 = 0; i2 < 7; i2++) {
        i3 = (i2 + 7 * i1) + 14 * i;
        in1_data[i3] = b_in1_data[i3];
      }
    }
  }

  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

void discretize_error_dynamics_FOH_RKV65_3DoF(const emlrtStack *sp, real_T N,
  const real_T tspan[2], const real_T x_ref[105], const real_T u_ref[30], real_T
  s_ref, real_T N_sub, real_T L, real_T b_I, real_T alpha, emxArray_real_T *A_k,
  emxArray_real_T *B_k_plus, emxArray_real_T *B_k_minus, emxArray_real_T *S_k,
  emxArray_real_T *d_k, emxArray_real_T *Delta)
{
  static const real_T dv2[36] = { 0.06, 0.019239962962962962,
    0.07669337037037037, 0.035975, 0.0, 0.107925, 1.3186834152331484, 0.0,
    -5.0420580636285619, 4.2206746483954136, -41.872591664327516, 0.0,
    159.4325621631375, -122.11921356501003, 5.531743066200054,
    -54.430156935316504, 0.0, 207.06725136501848, -158.61081378459,
    6.9918165859502421, -0.018597231062203234, -54.663741787281978, 0.0,
    207.95280625538936, -159.28895747449951, 7.0187437407969444,
    -0.018338785905045722, -0.00051194849978820987, 0.034389578683570357, 0.0,
    0.0, 0.25826245556335031, 0.4209371189673537, 4.40539646966931,
    -176.48311902429865, 172.36413340141507 };

  static const real_T dv[9] = { 0.0, 0.06, 0.095933333333333329, 0.1439, 0.4973,
    0.9725, 0.9995, 1.0, 1.0 };

  static const real_T dv1[9] = { 0.034389578683570357, 0.0, 0.0,
    0.25826245556335031, 0.4209371189673537, 4.40539646966931,
    -176.48311902429865, 172.36413340141507, 0.0 };

  __m128d r12;
  __m128d r14;
  emlrtStack b_st;
  emlrtStack st;
  emxArray_real_T r23;
  emxArray_real_T r24;
  emxArray_real_T *A_k_est;
  emxArray_real_T *A_kdot_ki;
  emxArray_real_T *B_k_minus_est;
  emxArray_real_T *B_k_minusdot_ki;
  emxArray_real_T *B_k_plus_est;
  emxArray_real_T *B_k_plusdot_ki;
  emxArray_real_T *S_k_est;
  emxArray_real_T *S_kdot_ki;
  emxArray_real_T *b_t_k;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  emxArray_real_T *r10;
  emxArray_real_T *r11;
  emxArray_real_T *r18;
  emxArray_real_T *r2;
  emxArray_real_T *r20;
  emxArray_real_T *r25;
  emxArray_real_T *r3;
  emxArray_real_T *r4;
  emxArray_real_T *r5;
  emxArray_real_T *r6;
  emxArray_real_T *r7;
  emxArray_real_T *r8;
  emxArray_real_T *r9;
  emxArray_real_T *t;
  emxArray_real_T *t_k;
  emxArray_real_T *x_est;
  emxArray_real_T *x_k;
  emxArray_real_T *xdot_ki;
  real_T u_ref_data[60];
  real_T d;
  real_T dt_sub;
  real_T dt_tmp;
  real_T *A_k_data;
  real_T *A_kdot_ki_data;
  real_T *B_k_minus_data;
  real_T *B_k_minusdot_ki_data;
  real_T *B_k_plus_data;
  real_T *B_k_plusdot_ki_data;
  real_T *S_k_data;
  real_T *S_kdot_ki_data;
  real_T *r13;
  real_T *r15;
  real_T *r16;
  real_T *r17;
  real_T *r19;
  real_T *r21;
  real_T *r22;
  real_T *t_data;
  real_T *t_k_data;
  real_T *x_est_data;
  real_T *x_k_data;
  real_T *xdot_ki_data;
  int32_T iv1[3];
  int32_T iv2[3];
  int32_T iv3[3];
  int32_T iv4[3];
  int32_T iv5[3];
  int32_T iv[2];
  int32_T b_i;
  int32_T b_loop_ub;
  int32_T b_loop_ub_tmp;
  int32_T c_loop_ub = 0;
  int32_T c_loop_ub_tmp;
  int32_T d_loop_ub_tmp;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  int32_T i5;
  int32_T i6;
  int32_T i7 = 0;
  int32_T i8;
  int32_T i9;
  int32_T j;
  int32_T k;
  int32_T loop_ub;
  int32_T loop_ub_tmp;
  int32_T n;
  int32_T nx;
  int32_T nx_tmp;
  int32_T scalarLB;
  int32_T vectorUB;
  (void)s_ref;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emlrtMEXProfilingFunctionEntry((char_T *)d_discretize_error_dynamics_FOH,
    isMexOutdated);

  /*  Discretization of the devition from a reference trajectory for a */
  /*  dynamical system assuming FOH control. Use Efficient Verner RK6(5) */
  /*  https://www.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT */
  /*  Should have local truncation error of  */
  /*  so N_sub = 1? Hopefully... */
  emlrtMEXProfilingStatement(4, isMexOutdated);
  emxInit_real_T(sp, &t, 2, &nc_emlrtRTEI);
  st.site = &emlrtRSI;
  linspace(&st, tspan[0], tspan[1], N, t);
  t_data = t->data;
  emlrtMEXProfilingStatement(5, isMexOutdated);
  if (t->size[1] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, t->size[1], &k_emlrtBCI, (emlrtConstCTX)
      sp);
  }

  if (t->size[1] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, t->size[1], &j_emlrtBCI, (emlrtConstCTX)
      sp);
  }

  dt_tmp = t_data[1] - t_data[0];
  emlrtMEXProfilingStatement(6, isMexOutdated);
  if (!(N - 1.0 >= 0.0)) {
    emlrtNonNegativeCheckR2012b(N - 1.0, &h_emlrtDCI, (emlrtConstCTX)sp);
  }

  i = (int32_T)muDoubleScalarFloor(N - 1.0);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &i_emlrtDCI, (emlrtConstCTX)sp);
  }

  i1 = A_k->size[0] * A_k->size[1] * A_k->size[2];
  A_k->size[0] = 7;
  A_k->size[1] = 7;
  loop_ub_tmp = (int32_T)(N - 1.0);
  A_k->size[2] = (int32_T)(N - 1.0);
  emxEnsureCapacity_real_T(sp, A_k, i1, &l_emlrtRTEI);
  A_k_data = A_k->data;
  b_loop_ub_tmp = 49 * (int32_T)(N - 1.0);
  for (i1 = 0; i1 < b_loop_ub_tmp; i1++) {
    A_k_data[i1] = 0.0;
  }

  emlrtMEXProfilingStatement(7, isMexOutdated);
  emlrtForLoopVectorCheckR2021a(1.0, 1.0, N - 1.0, mxDOUBLE_CLASS, (int32_T)(N -
    1.0), &e_emlrtRTEI, (emlrtConstCTX)sp);
  for (k = 0; k < loop_ub_tmp; k++) {
    emlrtMEXProfilingStatement(8, isMexOutdated);
    if (((int32_T)((uint32_T)k + 1U) < 1) || ((int32_T)((uint32_T)k + 1U) >
         A_k->size[2])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)k + 1U), 1, A_k->size[2],
        &l_emlrtBCI, (emlrtConstCTX)sp);
    }

    eye(&A_k_data[49 * k]);
    emlrtMEXProfilingStatement(9, isMexOutdated);
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }

  emlrtMEXProfilingStatement(10, isMexOutdated);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &j_emlrtDCI, (emlrtConstCTX)sp);
  }

  i1 = B_k_plus->size[0] * B_k_plus->size[1] * B_k_plus->size[2];
  B_k_plus->size[0] = 7;
  B_k_plus->size[1] = 2;
  B_k_plus->size[2] = (int32_T)(N - 1.0);
  emxEnsureCapacity_real_T(sp, B_k_plus, i1, &m_emlrtRTEI);
  B_k_plus_data = B_k_plus->data;
  b_loop_ub_tmp = 14 * (int32_T)(N - 1.0);
  for (i1 = 0; i1 < b_loop_ub_tmp; i1++) {
    B_k_plus_data[i1] = 0.0;
  }

  emlrtMEXProfilingStatement(11, isMexOutdated);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &k_emlrtDCI, (emlrtConstCTX)sp);
  }

  i1 = B_k_minus->size[0] * B_k_minus->size[1] * B_k_minus->size[2];
  B_k_minus->size[0] = 7;
  B_k_minus->size[1] = 2;
  B_k_minus->size[2] = (int32_T)(N - 1.0);
  emxEnsureCapacity_real_T(sp, B_k_minus, i1, &n_emlrtRTEI);
  B_k_minus_data = B_k_minus->data;
  for (i1 = 0; i1 < b_loop_ub_tmp; i1++) {
    B_k_minus_data[i1] = 0.0;
  }

  emlrtMEXProfilingStatement(12, isMexOutdated);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &l_emlrtDCI, (emlrtConstCTX)sp);
  }

  i1 = S_k->size[0] * S_k->size[1] * S_k->size[2];
  S_k->size[0] = 7;
  S_k->size[1] = 1;
  S_k->size[2] = (int32_T)(N - 1.0);
  emxEnsureCapacity_real_T(sp, S_k, i1, &o_emlrtRTEI);
  S_k_data = S_k->data;
  b_loop_ub_tmp = 7 * (int32_T)(N - 1.0);
  for (i1 = 0; i1 < b_loop_ub_tmp; i1++) {
    S_k_data[i1] = 0.0;
  }

  emlrtMEXProfilingStatement(13, isMexOutdated);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &m_emlrtDCI, (emlrtConstCTX)sp);
  }

  /* Delta = zeros([1, N - 1]); */
  /*     %% RKV65 */
  emlrtMEXProfilingStatement(14, isMexOutdated);
  if (t->size[1] - 1 < 1) {
    loop_ub = 0;
  } else {
    if (t->size[1] < 1) {
      emlrtDynamicBoundsCheckR2012b(1, 1, t->size[1], &i_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    if ((t->size[1] - 1 < 1) || (t->size[1] - 1 > t->size[1])) {
      emlrtDynamicBoundsCheckR2012b(t->size[1] - 1, 1, t->size[1], &h_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    loop_ub = t->size[1] - 1;
  }

  emlrtMEXProfilingStatement(15, isMexOutdated);
  if (t->size[1] < 2) {
    i1 = 0;
    i2 = 0;
  } else {
    i1 = 1;
    i2 = t->size[1];
  }

  emlrtMEXProfilingStatement(16, isMexOutdated);
  if (t->size[1] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, t->size[1], &g_emlrtBCI, (emlrtConstCTX)
      sp);
  }

  if (t->size[1] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, t->size[1], &f_emlrtBCI, (emlrtConstCTX)
      sp);
  }

  dt_sub = dt_tmp / N_sub;
  emlrtMEXProfilingStatement(17, isMexOutdated);
  if (N - 1.0 < 1.0) {
    b_loop_ub = 0;
  } else {
    if (N - 1.0 != i) {
      emlrtIntegerCheckR2012b(N - 1.0, &g_emlrtDCI, (emlrtConstCTX)sp);
    }

    if (((int32_T)(N - 1.0) < 1) || ((int32_T)(N - 1.0) > 15)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(N - 1.0), 1, 15, &e_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    b_loop_ub = (int32_T)(N - 1.0);
  }

  emxInit_real_T(sp, &x_k, 2, &p_emlrtRTEI);
  i3 = x_k->size[0] * x_k->size[1];
  x_k->size[0] = 7;
  x_k->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(sp, x_k, i3, &p_emlrtRTEI);
  x_k_data = x_k->data;
  for (i3 = 0; i3 < b_loop_ub; i3++) {
    for (i4 = 0; i4 < 7; i4++) {
      i5 = i4 + 7 * i3;
      x_k_data[i5] = x_ref[i5];
    }
  }

  emlrtMEXProfilingStatement(18, isMexOutdated);
  if (N - 1.0 < 1.0) {
    i3 = 0;
  } else {
    if (N - 1.0 != i) {
      emlrtIntegerCheckR2012b(N - 1.0, &f_emlrtDCI, (emlrtConstCTX)sp);
    }

    if (((int32_T)(N - 1.0) < 1) || ((int32_T)(N - 1.0) > 15)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(N - 1.0), 1, 15, &d_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    i3 = (int32_T)(N - 1.0);
  }

  emlrtMEXProfilingStatement(19, isMexOutdated);
  if (N < 2.0) {
    i4 = 0;
    i5 = 0;
  } else {
    i4 = 1;
    if (N != (int32_T)muDoubleScalarFloor(N)) {
      emlrtIntegerCheckR2012b(N, &e_emlrtDCI, (emlrtConstCTX)sp);
    }

    if (((int32_T)N < 1) || ((int32_T)N > 15)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)N, 1, 15, &c_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    i5 = (int32_T)N;
  }

  emlrtMEXProfilingStatement(20, isMexOutdated);
  emlrtMEXProfilingStatement(21, isMexOutdated);
  emlrtMEXProfilingStatement(22, isMexOutdated);

  /*  Propagation stage */
  emlrtMEXProfilingStatement(24, isMexOutdated);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &n_emlrtDCI, (emlrtConstCTX)sp);
  }

  emxInit_real_T(sp, &xdot_ki, 3, &q_emlrtRTEI);
  i6 = xdot_ki->size[0] * xdot_ki->size[1] * xdot_ki->size[2];
  xdot_ki->size[0] = 7;
  xdot_ki->size[1] = (int32_T)(N - 1.0);
  xdot_ki->size[2] = 9;
  emxEnsureCapacity_real_T(sp, xdot_ki, i6, &q_emlrtRTEI);
  xdot_ki_data = xdot_ki->data;
  c_loop_ub_tmp = b_loop_ub_tmp * 9;
  for (i6 = 0; i6 < c_loop_ub_tmp; i6++) {
    xdot_ki_data[i6] = 0.0;
  }

  emlrtMEXProfilingStatement(25, isMexOutdated);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &o_emlrtDCI, (emlrtConstCTX)sp);
  }

  emxInit_real_T(sp, &A_kdot_ki, 4, &r_emlrtRTEI);
  i6 = A_kdot_ki->size[0] * A_kdot_ki->size[1] * A_kdot_ki->size[2] *
    A_kdot_ki->size[3];
  A_kdot_ki->size[0] = 7;
  A_kdot_ki->size[1] = 7;
  c_loop_ub_tmp = (int8_T)(N - 1.0);
  A_kdot_ki->size[2] = (int8_T)(N - 1.0);
  A_kdot_ki->size[3] = 9;
  emxEnsureCapacity_real_T(sp, A_kdot_ki, i6, &r_emlrtRTEI);
  A_kdot_ki_data = A_kdot_ki->data;
  d_loop_ub_tmp = 49 * (int8_T)(N - 1.0) * 9;
  for (i6 = 0; i6 < d_loop_ub_tmp; i6++) {
    A_kdot_ki_data[i6] = 0.0;
  }

  emlrtMEXProfilingStatement(26, isMexOutdated);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &p_emlrtDCI, (emlrtConstCTX)sp);
  }

  emxInit_real_T(sp, &B_k_plusdot_ki, 4, &s_emlrtRTEI);
  i6 = B_k_plusdot_ki->size[0] * B_k_plusdot_ki->size[1] * B_k_plusdot_ki->size
    [2] * B_k_plusdot_ki->size[3];
  B_k_plusdot_ki->size[0] = 7;
  B_k_plusdot_ki->size[1] = 2;
  B_k_plusdot_ki->size[2] = (int8_T)(N - 1.0);
  B_k_plusdot_ki->size[3] = 9;
  emxEnsureCapacity_real_T(sp, B_k_plusdot_ki, i6, &s_emlrtRTEI);
  B_k_plusdot_ki_data = B_k_plusdot_ki->data;
  d_loop_ub_tmp = 14 * (int8_T)(N - 1.0) * 9;
  for (i6 = 0; i6 < d_loop_ub_tmp; i6++) {
    B_k_plusdot_ki_data[i6] = 0.0;
  }

  emlrtMEXProfilingStatement(27, isMexOutdated);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &q_emlrtDCI, (emlrtConstCTX)sp);
  }

  emxInit_real_T(sp, &B_k_minusdot_ki, 4, &t_emlrtRTEI);
  i6 = B_k_minusdot_ki->size[0] * B_k_minusdot_ki->size[1] *
    B_k_minusdot_ki->size[2] * B_k_minusdot_ki->size[3];
  B_k_minusdot_ki->size[0] = 7;
  B_k_minusdot_ki->size[1] = 2;
  B_k_minusdot_ki->size[2] = (int8_T)(N - 1.0);
  B_k_minusdot_ki->size[3] = 9;
  emxEnsureCapacity_real_T(sp, B_k_minusdot_ki, i6, &t_emlrtRTEI);
  B_k_minusdot_ki_data = B_k_minusdot_ki->data;
  for (i6 = 0; i6 < d_loop_ub_tmp; i6++) {
    B_k_minusdot_ki_data[i6] = 0.0;
  }

  emlrtMEXProfilingStatement(28, isMexOutdated);
  if (N - 1.0 != i) {
    emlrtIntegerCheckR2012b(N - 1.0, &r_emlrtDCI, (emlrtConstCTX)sp);
  }

  emxInit_real_T(sp, &S_kdot_ki, 4, &u_emlrtRTEI);
  i6 = S_kdot_ki->size[0] * S_kdot_ki->size[1] * S_kdot_ki->size[2] *
    S_kdot_ki->size[3];
  S_kdot_ki->size[0] = 7;
  S_kdot_ki->size[1] = 1;
  S_kdot_ki->size[2] = (int8_T)(N - 1.0);
  S_kdot_ki->size[3] = 9;
  emxEnsureCapacity_real_T(sp, S_kdot_ki, i6, &u_emlrtRTEI);
  S_kdot_ki_data = S_kdot_ki->data;
  d_loop_ub_tmp = 7 * (int8_T)(N - 1.0) * 9;
  for (i6 = 0; i6 < d_loop_ub_tmp; i6++) {
    S_kdot_ki_data[i6] = 0.0;
  }

  emlrtMEXProfilingStatement(29, isMexOutdated);
  emlrtMEXProfilingStatement(30, isMexOutdated);
  i6 = (int32_T)N_sub;
  emlrtForLoopVectorCheckR2021a(1.0, 1.0, N_sub, mxDOUBLE_CLASS, (int32_T)N_sub,
    &d_emlrtRTEI, (emlrtConstCTX)sp);
  emxInit_real_T(sp, &x_est, 2, &oc_emlrtRTEI);
  emxInit_real_T(sp, &A_k_est, 3, &pc_emlrtRTEI);
  emxInit_real_T(sp, &B_k_plus_est, 3, &qc_emlrtRTEI);
  emxInit_real_T(sp, &B_k_minus_est, 3, &rc_emlrtRTEI);
  emxInit_real_T(sp, &S_k_est, 3, &sc_emlrtRTEI);
  emxInit_real_T(sp, &t_k, 2, &v_emlrtRTEI);
  emxInit_real_T(sp, &r, 2, &cc_emlrtRTEI);
  emxInit_real_T(sp, &r1, 3, &fc_emlrtRTEI);
  emxInit_real_T(sp, &r2, 3, &hc_emlrtRTEI);
  emxInit_real_T(sp, &r3, 3, &tc_emlrtRTEI);
  emxInit_real_T(sp, &r4, 3, &tc_emlrtRTEI);
  emxInit_real_T(sp, &r5, 2, &bb_emlrtRTEI);
  emxInit_real_T(sp, &r6, 2, &bb_emlrtRTEI);
  emxInit_real_T(sp, &r7, 3, &lc_emlrtRTEI);
  emxInit_real_T(sp, &r8, 2, &w_emlrtRTEI);
  emxInit_real_T(sp, &r9, 2, &y_emlrtRTEI);
  emxInit_real_T(sp, &r10, 2, &x_emlrtRTEI);
  emxInit_real_T(sp, &b_t_k, 2, &db_emlrtRTEI);
  emxInit_real_T(sp, &r11, 2, &bb_emlrtRTEI);
  if (i6 - 1 >= 0) {
    scalarLB = (loop_ub / 2) << 1;
    vectorUB = scalarLB - 2;
    d = 0.0 * dt_sub;
    c_loop_ub = i2 - i1;
    i7 = i5 - i4;
  }

  for (k = 0; k < i6; k++) {
    real_T b_dt_sub;
    int32_T b_vectorUB;
    int32_T d_loop_ub;
    int32_T e_loop_ub_tmp;
    int32_T f_loop_ub_tmp;
    int32_T g_loop_ub_tmp;
    int32_T h_loop_ub_tmp;
    int32_T scalarLB_tmp;
    int32_T vectorUB_tmp;
    emlrtMEXProfilingStatement(31, isMexOutdated);
    emlrtMEXProfilingStatement(32, isMexOutdated);
    emlrtMEXProfilingStatement(33, isMexOutdated);
    emlrtMEXProfilingStatement(34, isMexOutdated);
    emlrtMEXProfilingStatement(35, isMexOutdated);
    emlrtMEXProfilingStatement(36, isMexOutdated);
    b_dt_sub = dt_sub * (((real_T)k + 1.0) - 1.0);
    i8 = t_k->size[0] * t_k->size[1];
    t_k->size[0] = 1;
    t_k->size[1] = loop_ub;
    emxEnsureCapacity_real_T(sp, t_k, i8, &v_emlrtRTEI);
    t_k_data = t_k->data;
    for (i8 = 0; i8 <= vectorUB; i8 += 2) {
      r12 = _mm_loadu_pd(&t_data[i8]);
      _mm_storeu_pd(&t_k_data[i8], _mm_add_pd(r12, _mm_set1_pd(b_dt_sub)));
    }

    for (i8 = scalarLB; i8 < loop_ub; i8++) {
      t_k_data[i8] = t_data[i8] + b_dt_sub;
    }

    emlrtMEXProfilingStatement(37, isMexOutdated);
    i8 = r8->size[0] * r8->size[1];
    r8->size[0] = 1;
    b_loop_ub = t_k->size[1];
    r8->size[1] = t_k->size[1];
    emxEnsureCapacity_real_T(sp, r8, i8, &w_emlrtRTEI);
    r13 = r8->data;
    n = (t_k->size[1] / 2) << 1;
    b_vectorUB = n - 2;
    for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
      r12 = _mm_loadu_pd(&t_k_data[i8]);
      r14 = _mm_loadu_pd(&t_data[i8]);
      _mm_storeu_pd(&r13[i8], _mm_div_pd(_mm_sub_pd(_mm_add_pd(r12, _mm_set1_pd
        (d)), r14), _mm_set1_pd(dt_tmp)));
    }

    for (i8 = n; i8 < b_loop_ub; i8++) {
      r13[i8] = ((t_k_data[i8] + d) - t_data[i8]) / dt_tmp;
    }

    st.site = &b_emlrtRSI;
    nx_tmp = r8->size[1];
    n = 1;
    if (r8->size[1] > 1) {
      n = r8->size[1];
    }

    if (r8->size[1] > muIntScalarMax_sint32(nx_tmp, n)) {
      emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
        "Coder:toolbox:reshape_emptyReshapeLimit",
        "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }

    emlrtMEXProfilingStatement(38, isMexOutdated);
    if ((c_loop_ub != t_k->size[1]) && ((c_loop_ub != 1) && (t_k->size[1] != 1)))
    {
      emlrtDimSizeImpxCheckR2021b(c_loop_ub, t_k->size[1], &ib_emlrtECI,
        (emlrtConstCTX)sp);
    }

    if (c_loop_ub == t_k->size[1]) {
      i8 = r10->size[0] * r10->size[1];
      r10->size[0] = 1;
      r10->size[1] = c_loop_ub;
      emxEnsureCapacity_real_T(sp, r10, i8, &x_emlrtRTEI);
      r15 = r10->data;
      n = (c_loop_ub / 2) << 1;
      b_vectorUB = n - 2;
      for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
        r12 = _mm_loadu_pd(&t_data[i1 + i8]);
        r14 = _mm_loadu_pd(&t_k_data[i8]);
        _mm_storeu_pd(&r15[i8], _mm_sub_pd(r12, r14));
      }

      for (i8 = n; i8 < c_loop_ub; i8++) {
        r15[i8] = t_data[i1 + i8] - t_k_data[i8];
      }
    } else {
      st.site = &c_emlrtRSI;
      binary_expand_op_6(&st, r10, t, i1, i2, t_k);
      r15 = r10->data;
    }

    i8 = r9->size[0] * r9->size[1];
    r9->size[0] = 1;
    b_loop_ub = r10->size[1];
    r9->size[1] = r10->size[1];
    emxEnsureCapacity_real_T(sp, r9, i8, &y_emlrtRTEI);
    r16 = r9->data;
    scalarLB_tmp = (r10->size[1] / 2) << 1;
    vectorUB_tmp = scalarLB_tmp - 2;
    for (i8 = 0; i8 <= vectorUB_tmp; i8 += 2) {
      r12 = _mm_loadu_pd(&r15[i8]);
      _mm_storeu_pd(&r16[i8], _mm_div_pd(_mm_sub_pd(r12, _mm_set1_pd(d)),
        _mm_set1_pd(dt_tmp)));
    }

    for (i8 = scalarLB_tmp; i8 < b_loop_ub; i8++) {
      r16[i8] = (r15[i8] - d) / dt_tmp;
    }

    st.site = &c_emlrtRSI;
    nx = r9->size[1];
    n = 1;
    if (r9->size[1] > 1) {
      n = r9->size[1];
    }

    if (r9->size[1] > muIntScalarMax_sint32(nx, n)) {
      emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
        "Coder:toolbox:reshape_emptyReshapeLimit",
        "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }

    emlrtMEXProfilingStatement(39, isMexOutdated);
    emlrtMEXProfilingStatement(40, isMexOutdated);
    if ((r8->size[1] != i3) && ((r8->size[1] != 1) && (i3 != 1))) {
      emlrtDimSizeImpxCheckR2021b(r8->size[1], i3, &fb_emlrtECI, (emlrtConstCTX)
        sp);
    }

    d_loop_ub = t_k->size[1];
    if ((c_loop_ub != t_k->size[1]) && ((c_loop_ub != 1) && (t_k->size[1] != 1)))
    {
      emlrtDimSizeImpxCheckR2021b(c_loop_ub, t_k->size[1], &hb_emlrtECI,
        (emlrtConstCTX)sp);
    }

    if ((r9->size[1] != i7) && ((r9->size[1] != 1) && (i7 != 1))) {
      emlrtDimSizeImpxCheckR2021b(r9->size[1], i7, &gb_emlrtECI, (emlrtConstCTX)
        sp);
    }

    if (r8->size[1] == i3) {
      i8 = r5->size[0] * r5->size[1];
      r5->size[0] = 2;
      r5->size[1] = r8->size[1];
      emxEnsureCapacity_real_T(sp, r5, i8, &bb_emlrtRTEI);
      r17 = r5->data;
      for (i8 = 0; i8 < nx_tmp; i8++) {
        _mm_storeu_pd(&r17[2 * i8], _mm_mul_pd(_mm_set1_pd(r13[i8]),
          _mm_loadu_pd(&u_ref[i8 << 1])));
      }
    } else {
      st.site = &ub_emlrtRSI;
      binary_expand_op_2(&st, r5, r8, u_ref, i3);
      r17 = r5->data;
    }

    if (r9->size[1] == i7) {
      i8 = r6->size[0] * r6->size[1];
      r6->size[0] = 2;
      r6->size[1] = r9->size[1];
      emxEnsureCapacity_real_T(sp, r6, i8, &cb_emlrtRTEI);
      r13 = r6->data;
      for (i8 = 0; i8 < nx; i8++) {
        _mm_storeu_pd(&r13[2 * i8], _mm_mul_pd(_mm_set1_pd(r16[i8]),
          _mm_loadu_pd(&u_ref[(i4 + i8) << 1])));
      }
    } else {
      st.site = &wb_emlrtRSI;
      binary_expand_op_1(&st, r6, r9, u_ref, i4, i5);
      r13 = r6->data;
    }

    if ((r5->size[1] != r6->size[1]) && ((r5->size[1] != 1) && (r6->size[1] != 1)))
    {
      emlrtDimSizeImpxCheckR2021b(r5->size[1], r6->size[1], &fb_emlrtECI,
        (emlrtConstCTX)sp);
    }

    emlrtMEXProfilingStatement(41, isMexOutdated);
    if (r5->size[1] == r6->size[1]) {
      iv1[2] = r8->size[1];
      nx = r9->size[1];
      i8 = b_t_k->size[0] * b_t_k->size[1];
      b_t_k->size[0] = 1;
      b_t_k->size[1] = t_k->size[1];
      emxEnsureCapacity_real_T(sp, b_t_k, i8, &db_emlrtRTEI);
      x_est_data = b_t_k->data;
      n = (t_k->size[1] / 2) << 1;
      b_vectorUB = n - 2;
      for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
        r12 = _mm_loadu_pd(&t_k_data[i8]);
        _mm_storeu_pd(&x_est_data[i8], _mm_add_pd(r12, _mm_set1_pd(d)));
      }

      for (i8 = n; i8 < d_loop_ub; i8++) {
        x_est_data[i8] = t_k_data[i8] + d;
      }

      i8 = r11->size[0] * r11->size[1];
      r11->size[0] = 2;
      r11->size[1] = r5->size[1];
      emxEnsureCapacity_real_T(sp, r11, i8, &bb_emlrtRTEI);
      r16 = r11->data;
      d_loop_ub_tmp = r5->size[1] << 1;
      n = (d_loop_ub_tmp / 2) << 1;
      b_vectorUB = n - 2;
      for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
        r12 = _mm_loadu_pd(&r17[i8]);
        r14 = _mm_loadu_pd(&r13[i8]);
        _mm_storeu_pd(&r16[i8], _mm_add_pd(r12, r14));
      }

      for (i8 = n; i8 < d_loop_ub_tmp; i8++) {
        r16[i8] = r17[i8] + r13[i8];
      }

      r23 = *r8;
      iv2[0] = 1;
      iv2[1] = 1;
      iv2[2] = iv1[2];
      r23.size = &iv2[0];
      r23.numDimensions = 3;
      r24 = *r9;
      iv3[0] = 1;
      iv3[1] = 1;
      iv3[2] = nx;
      r24.size = &iv3[0];
      r24.numDimensions = 3;
      st.site = &d_emlrtRSI;
      STM_diff_eq_FOH(&st, b_t_k, x_k, r11, &r23, &r24, A_k, B_k_plus, B_k_minus,
                      S_k, L, b_I, alpha, r, r1, r2, r3, r4);
      r13 = r4->data;
      r16 = r3->data;
      r19 = r2->data;
      r21 = r1->data;
      r22 = r->data;
    } else {
      st.site = &d_emlrtRSI;
      binary_expand_op_3(&st, d_emlrtRSI, t_k, dt_sub, x_k, r5, r6, r8, r9, A_k,
                         B_k_plus, B_k_minus, S_k, L, b_I, alpha, r, r1, r2, r3,
                         r4);
      r13 = r4->data;
      r16 = r3->data;
      r19 = r2->data;
      r21 = r1->data;
      r22 = r->data;
    }

    iv[0] = 7;
    iv[1] = (int32_T)(N - 1.0);
    emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &r->size[0], 2, &ab_emlrtECI,
      (emlrtCTX)sp);
    for (i8 = 0; i8 < loop_ub_tmp; i8++) {
      for (nx_tmp = 0; nx_tmp < 7; nx_tmp++) {
        i9 = nx_tmp + 7 * i8;
        xdot_ki_data[i9] = r22[i9];
      }
    }

    iv1[0] = 7;
    iv1[1] = 7;
    iv1[2] = (int8_T)(N - 1.0);
    emlrtSubAssignSizeCheckR2012b(&iv1[0], 3, &r1->size[0], 3, &x_emlrtECI,
      (emlrtCTX)sp);
    for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
      for (nx_tmp = 0; nx_tmp < 7; nx_tmp++) {
        for (i9 = 0; i9 < 7; i9++) {
          n = (i9 + 7 * nx_tmp) + 49 * i8;
          A_kdot_ki_data[n] = r21[n];
        }
      }
    }

    iv1[0] = 7;
    iv1[1] = 2;
    iv1[2] = (int8_T)(N - 1.0);
    emlrtSubAssignSizeCheckR2012b(&iv1[0], 3, &r2->size[0], 3, &v_emlrtECI,
      (emlrtCTX)sp);
    for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
      for (nx_tmp = 0; nx_tmp < 2; nx_tmp++) {
        for (i9 = 0; i9 < 7; i9++) {
          n = (i9 + 7 * nx_tmp) + 14 * i8;
          B_k_plusdot_ki_data[n] = r19[n];
        }
      }
    }

    iv1[0] = 7;
    iv1[1] = 2;
    iv1[2] = (int8_T)(N - 1.0);
    emlrtSubAssignSizeCheckR2012b(&iv1[0], 3, &r3->size[0], 3, &t_emlrtECI,
      (emlrtCTX)sp);
    for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
      for (nx_tmp = 0; nx_tmp < 2; nx_tmp++) {
        for (i9 = 0; i9 < 7; i9++) {
          n = (i9 + 7 * nx_tmp) + 14 * i8;
          B_k_minusdot_ki_data[n] = r16[n];
        }
      }
    }

    iv1[0] = 7;
    iv1[1] = 1;
    iv1[2] = (int8_T)(N - 1.0);
    emlrtSubAssignSizeCheckR2012b(&iv1[0], 3, &r4->size[0], 3, &r_emlrtECI,
      (emlrtCTX)sp);
    for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
      for (nx_tmp = 0; nx_tmp < 7; nx_tmp++) {
        i9 = nx_tmp + 7 * i8;
        S_kdot_ki_data[i9] = r13[i9];
      }
    }

    emlrtMEXProfilingStatement(42, isMexOutdated);
    d_loop_ub_tmp = 7 * x_k->size[1];
    e_loop_ub_tmp = 49 * A_k->size[2];
    f_loop_ub_tmp = 14 * B_k_plus->size[2];
    g_loop_ub_tmp = 14 * B_k_minus->size[2];
    h_loop_ub_tmp = 7 * S_k->size[2];
    for (b_i = 0; b_i < 8; b_i++) {
      int32_T b_scalarLB_tmp;
      int32_T b_vectorUB_tmp;
      int32_T i_loop_ub_tmp;
      emlrtMEXProfilingStatement(43, isMexOutdated);
      b_dt_sub = dv[b_i + 1] * dt_sub;
      i8 = r8->size[0] * r8->size[1];
      r8->size[0] = 1;
      i_loop_ub_tmp = t_k->size[1];
      r8->size[1] = t_k->size[1];
      emxEnsureCapacity_real_T(sp, r8, i8, &gb_emlrtRTEI);
      r13 = r8->data;
      b_scalarLB_tmp = (t_k->size[1] / 2) << 1;
      b_vectorUB_tmp = b_scalarLB_tmp - 2;
      for (i8 = 0; i8 <= b_vectorUB_tmp; i8 += 2) {
        r12 = _mm_loadu_pd(&t_k_data[i8]);
        r14 = _mm_loadu_pd(&t_data[i8]);
        _mm_storeu_pd(&r13[i8], _mm_div_pd(_mm_sub_pd(_mm_add_pd(r12,
          _mm_set1_pd(b_dt_sub)), r14), _mm_set1_pd(dt_tmp)));
      }

      for (i8 = b_scalarLB_tmp; i8 < i_loop_ub_tmp; i8++) {
        r13[i8] = ((t_k_data[i8] + b_dt_sub) - t_data[i8]) / dt_tmp;
      }

      st.site = &e_emlrtRSI;
      nx_tmp = r8->size[1];
      n = 1;
      if (r8->size[1] > 1) {
        n = r8->size[1];
      }

      if (r8->size[1] > muIntScalarMax_sint32(nx_tmp, n)) {
        emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
          "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
      }

      emlrtMEXProfilingStatement(44, isMexOutdated);
      if ((c_loop_ub != t_k->size[1]) && ((c_loop_ub != 1) && (t_k->size[1] != 1)))
      {
        emlrtDimSizeImpxCheckR2021b(c_loop_ub, t_k->size[1], &eb_emlrtECI,
          (emlrtConstCTX)sp);
      }

      i8 = r9->size[0] * r9->size[1];
      r9->size[0] = 1;
      r9->size[1] = b_loop_ub;
      emxEnsureCapacity_real_T(sp, r9, i8, &mb_emlrtRTEI);
      r16 = r9->data;
      for (i8 = 0; i8 <= vectorUB_tmp; i8 += 2) {
        r12 = _mm_loadu_pd(&r15[i8]);
        _mm_storeu_pd(&r16[i8], _mm_div_pd(_mm_sub_pd(r12, _mm_set1_pd(b_dt_sub)),
          _mm_set1_pd(dt_tmp)));
      }

      for (i8 = scalarLB_tmp; i8 < b_loop_ub; i8++) {
        r16[i8] = (r15[i8] - b_dt_sub) / dt_tmp;
      }

      st.site = &f_emlrtRSI;
      nx = r9->size[1];
      n = 1;
      if (r9->size[1] > 1) {
        n = r9->size[1];
      }

      if (r9->size[1] > muIntScalarMax_sint32(nx, n)) {
        emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
          "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
      }

      emlrtMEXProfilingStatement(45, isMexOutdated);
      emlrtMEXProfilingStatement(46, isMexOutdated);
      if ((r8->size[1] != i3) && ((r8->size[1] != 1) && (i3 != 1))) {
        emlrtDimSizeImpxCheckR2021b(r8->size[1], i3, &bb_emlrtECI,
          (emlrtConstCTX)sp);
      }

      if ((c_loop_ub != t_k->size[1]) && ((c_loop_ub != 1) && (t_k->size[1] != 1)))
      {
        emlrtDimSizeImpxCheckR2021b(c_loop_ub, t_k->size[1], &db_emlrtECI,
          (emlrtConstCTX)sp);
      }

      if ((r9->size[1] != i7) && ((r9->size[1] != 1) && (i7 != 1))) {
        emlrtDimSizeImpxCheckR2021b(r9->size[1], i7, &cb_emlrtECI,
          (emlrtConstCTX)sp);
      }

      if (r8->size[1] == i3) {
        i8 = r5->size[0] * r5->size[1];
        r5->size[0] = 2;
        r5->size[1] = r8->size[1];
        emxEnsureCapacity_real_T(sp, r5, i8, &rb_emlrtRTEI);
        r17 = r5->data;
        for (i8 = 0; i8 < nx_tmp; i8++) {
          _mm_storeu_pd(&r17[2 * i8], _mm_mul_pd(_mm_set1_pd(r13[i8]),
            _mm_loadu_pd(&u_ref[i8 << 1])));
        }
      } else {
        st.site = &tb_emlrtRSI;
        binary_expand_op_2(&st, r5, r8, u_ref, i3);
        r17 = r5->data;
      }

      if (r9->size[1] == i7) {
        i8 = r6->size[0] * r6->size[1];
        r6->size[0] = 2;
        r6->size[1] = r9->size[1];
        emxEnsureCapacity_real_T(sp, r6, i8, &sb_emlrtRTEI);
        r13 = r6->data;
        for (i8 = 0; i8 < nx; i8++) {
          _mm_storeu_pd(&r13[2 * i8], _mm_mul_pd(_mm_set1_pd(r16[i8]),
            _mm_loadu_pd(&u_ref[(i4 + i8) << 1])));
        }
      } else {
        st.site = &vb_emlrtRSI;
        binary_expand_op_1(&st, r6, r9, u_ref, i4, i5);
        r13 = r6->data;
      }

      if ((r5->size[1] != r6->size[1]) && ((r5->size[1] != 1) && (r6->size[1] !=
            1))) {
        emlrtDimSizeImpxCheckR2021b(r5->size[1], r6->size[1], &bb_emlrtECI,
          (emlrtConstCTX)sp);
      }

      emlrtMEXProfilingStatement(47, isMexOutdated);
      i8 = x_est->size[0] * x_est->size[1];
      x_est->size[0] = 7;
      x_est->size[1] = x_k->size[1];
      emxEnsureCapacity_real_T(sp, x_est, i8, &ub_emlrtRTEI);
      x_est_data = x_est->data;
      for (i8 = 0; i8 < d_loop_ub_tmp; i8++) {
        x_est_data[i8] = x_k_data[i8];
      }

      emlrtMEXProfilingStatement(48, isMexOutdated);
      i8 = A_k_est->size[0] * A_k_est->size[1] * A_k_est->size[2];
      A_k_est->size[0] = 7;
      A_k_est->size[1] = 7;
      A_k_est->size[2] = A_k->size[2];
      emxEnsureCapacity_real_T(sp, A_k_est, i8, &wb_emlrtRTEI);
      x_est_data = A_k_est->data;
      for (i8 = 0; i8 < e_loop_ub_tmp; i8++) {
        x_est_data[i8] = A_k_data[i8];
      }

      emlrtMEXProfilingStatement(49, isMexOutdated);
      i8 = B_k_plus_est->size[0] * B_k_plus_est->size[1] * B_k_plus_est->size[2];
      B_k_plus_est->size[0] = 7;
      B_k_plus_est->size[1] = 2;
      B_k_plus_est->size[2] = B_k_plus->size[2];
      emxEnsureCapacity_real_T(sp, B_k_plus_est, i8, &yb_emlrtRTEI);
      x_est_data = B_k_plus_est->data;
      for (i8 = 0; i8 < f_loop_ub_tmp; i8++) {
        x_est_data[i8] = B_k_plus_data[i8];
      }

      emlrtMEXProfilingStatement(50, isMexOutdated);
      i8 = B_k_minus_est->size[0] * B_k_minus_est->size[1] * B_k_minus_est->
        size[2];
      B_k_minus_est->size[0] = 7;
      B_k_minus_est->size[1] = 2;
      B_k_minus_est->size[2] = B_k_minus->size[2];
      emxEnsureCapacity_real_T(sp, B_k_minus_est, i8, &ac_emlrtRTEI);
      x_est_data = B_k_minus_est->data;
      for (i8 = 0; i8 < g_loop_ub_tmp; i8++) {
        x_est_data[i8] = B_k_minus_data[i8];
      }

      emlrtMEXProfilingStatement(51, isMexOutdated);
      i8 = S_k_est->size[0] * S_k_est->size[1] * S_k_est->size[2];
      S_k_est->size[0] = 7;
      S_k_est->size[1] = 1;
      S_k_est->size[2] = S_k->size[2];
      emxEnsureCapacity_real_T(sp, S_k_est, i8, &bc_emlrtRTEI);
      x_est_data = S_k_est->data;
      for (i8 = 0; i8 < h_loop_ub_tmp; i8++) {
        x_est_data[i8] = S_k_data[i8];
      }

      emlrtMEXProfilingStatement(52, isMexOutdated);
      for (j = 0; j <= b_i; j++) {
        real_T d1;
        real_T d2;
        emlrtMEXProfilingStatement(53, isMexOutdated);
        d1 = (real_T)(b_i * (b_i + 1)) / 2.0 + ((real_T)j + 1.0);
        i8 = (int32_T)muDoubleScalarFloor(d1);
        if (d1 != i8) {
          emlrtIntegerCheckR2012b(d1, &s_emlrtDCI, (emlrtConstCTX)sp);
        }

        nx_tmp = (int32_T)d1;
        d2 = dv2[nx_tmp - 1] * dt_sub;
        i9 = r->size[0] * r->size[1];
        r->size[0] = 7;
        r->size[1] = (int32_T)(N - 1.0);
        emxEnsureCapacity_real_T(sp, r, i9, &cc_emlrtRTEI);
        r22 = r->data;
        for (i9 = 0; i9 < loop_ub_tmp; i9++) {
          r12 = _mm_loadu_pd(&xdot_ki_data[7 * i9 + 7 * xdot_ki->size[1] * j]);
          r14 = _mm_set1_pd(d2);
          _mm_storeu_pd(&r22[7 * i9], _mm_mul_pd(r14, r12));
          r12 = _mm_loadu_pd(&xdot_ki_data[(7 * i9 + 7 * xdot_ki->size[1] * j) +
                             2]);
          _mm_storeu_pd(&r22[7 * i9 + 2], _mm_mul_pd(r14, r12));
          r12 = _mm_loadu_pd(&xdot_ki_data[(7 * i9 + 7 * xdot_ki->size[1] * j) +
                             4]);
          _mm_storeu_pd(&r22[7 * i9 + 4], _mm_mul_pd(r14, r12));
          r22[7 * i9 + 6] = d2 * xdot_ki_data[(7 * i9 + 7 * xdot_ki->size[1] * j)
            + 6];
        }

        if ((r->size[1] != x_est->size[1]) && ((r->size[1] != 1) && (x_est->
              size[1] != 1))) {
          emlrtDimSizeImpxCheckR2021b(r->size[1], x_est->size[1], &y_emlrtECI,
            (emlrtConstCTX)sp);
        }

        if (r->size[1] == x_est->size[1]) {
          d_loop_ub = 7 * r->size[1];
          i9 = x_est->size[0] * x_est->size[1];
          x_est->size[0] = 7;
          x_est->size[1] = r->size[1];
          emxEnsureCapacity_real_T(sp, x_est, i9, &ec_emlrtRTEI);
          x_est_data = x_est->data;
          n = (d_loop_ub / 2) << 1;
          b_vectorUB = n - 2;
          for (i9 = 0; i9 <= b_vectorUB; i9 += 2) {
            r12 = _mm_loadu_pd(&r22[i9]);
            r14 = _mm_loadu_pd(&x_est_data[i9]);
            _mm_storeu_pd(&x_est_data[i9], _mm_add_pd(r12, r14));
          }

          for (i9 = n; i9 < d_loop_ub; i9++) {
            x_est_data[i9] += r22[i9];
          }
        } else {
          st.site = &mb_emlrtRSI;
          d_plus(&st, x_est, r);
        }

        emlrtMEXProfilingStatement(54, isMexOutdated);
        if (nx_tmp != i8) {
          emlrtIntegerCheckR2012b(d1, &t_emlrtDCI, (emlrtConstCTX)sp);
        }

        i9 = r1->size[0] * r1->size[1] * r1->size[2];
        r1->size[0] = 7;
        r1->size[1] = 7;
        r1->size[2] = (int8_T)(N - 1.0);
        emxEnsureCapacity_real_T(sp, r1, i9, &fc_emlrtRTEI);
        r21 = r1->data;
        for (i9 = 0; i9 < c_loop_ub_tmp; i9++) {
          r12 = _mm_set1_pd(d2);
          for (n = 0; n < 7; n++) {
            nx = 7 * n + 49 * i9;
            r14 = _mm_loadu_pd(&A_kdot_ki_data[nx + 49 * A_kdot_ki->size[2] * j]);
            _mm_storeu_pd(&r21[nx], _mm_mul_pd(r12, r14));
            r14 = _mm_loadu_pd(&A_kdot_ki_data[((7 * n + 49 * i9) + 49 *
              A_kdot_ki->size[2] * j) + 2]);
            _mm_storeu_pd(&r21[nx + 2], _mm_mul_pd(r12, r14));
            r14 = _mm_loadu_pd(&A_kdot_ki_data[((7 * n + 49 * i9) + 49 *
              A_kdot_ki->size[2] * j) + 4]);
            _mm_storeu_pd(&r21[nx + 4], _mm_mul_pd(r12, r14));
            r21[nx + 6] = d2 * A_kdot_ki_data[((7 * n + 49 * i9) + 49 *
              A_kdot_ki->size[2] * j) + 6];
          }
        }

        if ((r1->size[2] != A_k_est->size[2]) && ((r1->size[2] != 1) &&
             (A_k_est->size[2] != 1))) {
          emlrtDimSizeImpxCheckR2021b(r1->size[2], A_k_est->size[2], &w_emlrtECI,
            (emlrtConstCTX)sp);
        }

        if (r1->size[2] == A_k_est->size[2]) {
          d_loop_ub = 49 * r1->size[2];
          i9 = A_k_est->size[0] * A_k_est->size[1] * A_k_est->size[2];
          A_k_est->size[0] = 7;
          A_k_est->size[1] = 7;
          A_k_est->size[2] = r1->size[2];
          emxEnsureCapacity_real_T(sp, A_k_est, i9, &gc_emlrtRTEI);
          x_est_data = A_k_est->data;
          n = (d_loop_ub / 2) << 1;
          b_vectorUB = n - 2;
          for (i9 = 0; i9 <= b_vectorUB; i9 += 2) {
            r12 = _mm_loadu_pd(&r21[i9]);
            r14 = _mm_loadu_pd(&x_est_data[i9]);
            _mm_storeu_pd(&x_est_data[i9], _mm_add_pd(r12, r14));
          }

          for (i9 = n; i9 < d_loop_ub; i9++) {
            x_est_data[i9] += r21[i9];
          }
        } else {
          st.site = &lb_emlrtRSI;
          c_plus(&st, A_k_est, r1);
        }

        emlrtMEXProfilingStatement(55, isMexOutdated);
        if (nx_tmp != i8) {
          emlrtIntegerCheckR2012b(d1, &u_emlrtDCI, (emlrtConstCTX)sp);
        }

        i9 = r2->size[0] * r2->size[1] * r2->size[2];
        r2->size[0] = 7;
        r2->size[1] = 2;
        r2->size[2] = (int8_T)(N - 1.0);
        emxEnsureCapacity_real_T(sp, r2, i9, &hc_emlrtRTEI);
        r19 = r2->data;
        for (i9 = 0; i9 < c_loop_ub_tmp; i9++) {
          r12 = _mm_set1_pd(d2);
          for (n = 0; n < 2; n++) {
            nx = 7 * n + 14 * i9;
            r14 = _mm_loadu_pd(&B_k_plusdot_ki_data[nx + 14 *
                               B_k_plusdot_ki->size[2] * j]);
            _mm_storeu_pd(&r19[nx], _mm_mul_pd(r12, r14));
            r14 = _mm_loadu_pd(&B_k_plusdot_ki_data[((7 * n + 14 * i9) + 14 *
              B_k_plusdot_ki->size[2] * j) + 2]);
            _mm_storeu_pd(&r19[nx + 2], _mm_mul_pd(r12, r14));
            r14 = _mm_loadu_pd(&B_k_plusdot_ki_data[((7 * n + 14 * i9) + 14 *
              B_k_plusdot_ki->size[2] * j) + 4]);
            _mm_storeu_pd(&r19[nx + 4], _mm_mul_pd(r12, r14));
            r19[nx + 6] = d2 * B_k_plusdot_ki_data[((7 * n + 14 * i9) + 14 *
              B_k_plusdot_ki->size[2] * j) + 6];
          }
        }

        if ((r2->size[2] != B_k_plus_est->size[2]) && ((r2->size[2] != 1) &&
             (B_k_plus_est->size[2] != 1))) {
          emlrtDimSizeImpxCheckR2021b(r2->size[2], B_k_plus_est->size[2],
            &u_emlrtECI, (emlrtConstCTX)sp);
        }

        if (r2->size[2] == B_k_plus_est->size[2]) {
          d_loop_ub = 14 * r2->size[2];
          i9 = B_k_plus_est->size[0] * B_k_plus_est->size[1] *
            B_k_plus_est->size[2];
          B_k_plus_est->size[0] = 7;
          B_k_plus_est->size[1] = 2;
          B_k_plus_est->size[2] = r2->size[2];
          emxEnsureCapacity_real_T(sp, B_k_plus_est, i9, &ic_emlrtRTEI);
          x_est_data = B_k_plus_est->data;
          n = (d_loop_ub / 2) << 1;
          b_vectorUB = n - 2;
          for (i9 = 0; i9 <= b_vectorUB; i9 += 2) {
            r12 = _mm_loadu_pd(&r19[i9]);
            r14 = _mm_loadu_pd(&x_est_data[i9]);
            _mm_storeu_pd(&x_est_data[i9], _mm_add_pd(r12, r14));
          }

          for (i9 = n; i9 < d_loop_ub; i9++) {
            x_est_data[i9] += r19[i9];
          }
        } else {
          st.site = &kb_emlrtRSI;
          b_plus(&st, B_k_plus_est, r2);
        }

        emlrtMEXProfilingStatement(56, isMexOutdated);
        if (nx_tmp != i8) {
          emlrtIntegerCheckR2012b(d1, &v_emlrtDCI, (emlrtConstCTX)sp);
        }

        i9 = r2->size[0] * r2->size[1] * r2->size[2];
        r2->size[0] = 7;
        r2->size[1] = 2;
        r2->size[2] = (int8_T)(N - 1.0);
        emxEnsureCapacity_real_T(sp, r2, i9, &jc_emlrtRTEI);
        r19 = r2->data;
        for (i9 = 0; i9 < c_loop_ub_tmp; i9++) {
          r12 = _mm_set1_pd(d2);
          for (n = 0; n < 2; n++) {
            nx = 7 * n + 14 * i9;
            r14 = _mm_loadu_pd(&B_k_minusdot_ki_data[nx + 14 *
                               B_k_minusdot_ki->size[2] * j]);
            _mm_storeu_pd(&r19[nx], _mm_mul_pd(r12, r14));
            r14 = _mm_loadu_pd(&B_k_minusdot_ki_data[((7 * n + 14 * i9) + 14 *
              B_k_minusdot_ki->size[2] * j) + 2]);
            _mm_storeu_pd(&r19[nx + 2], _mm_mul_pd(r12, r14));
            r14 = _mm_loadu_pd(&B_k_minusdot_ki_data[((7 * n + 14 * i9) + 14 *
              B_k_minusdot_ki->size[2] * j) + 4]);
            _mm_storeu_pd(&r19[nx + 4], _mm_mul_pd(r12, r14));
            r19[nx + 6] = d2 * B_k_minusdot_ki_data[((7 * n + 14 * i9) + 14 *
              B_k_minusdot_ki->size[2] * j) + 6];
          }
        }

        if ((r2->size[2] != B_k_minus_est->size[2]) && ((r2->size[2] != 1) &&
             (B_k_minus_est->size[2] != 1))) {
          emlrtDimSizeImpxCheckR2021b(r2->size[2], B_k_minus_est->size[2],
            &s_emlrtECI, (emlrtConstCTX)sp);
        }

        if (r2->size[2] == B_k_minus_est->size[2]) {
          d_loop_ub = 14 * r2->size[2];
          i9 = B_k_minus_est->size[0] * B_k_minus_est->size[1] *
            B_k_minus_est->size[2];
          B_k_minus_est->size[0] = 7;
          B_k_minus_est->size[1] = 2;
          B_k_minus_est->size[2] = r2->size[2];
          emxEnsureCapacity_real_T(sp, B_k_minus_est, i9, &kc_emlrtRTEI);
          x_est_data = B_k_minus_est->data;
          n = (d_loop_ub / 2) << 1;
          b_vectorUB = n - 2;
          for (i9 = 0; i9 <= b_vectorUB; i9 += 2) {
            r12 = _mm_loadu_pd(&r19[i9]);
            r14 = _mm_loadu_pd(&x_est_data[i9]);
            _mm_storeu_pd(&x_est_data[i9], _mm_add_pd(r12, r14));
          }

          for (i9 = n; i9 < d_loop_ub; i9++) {
            x_est_data[i9] += r19[i9];
          }
        } else {
          st.site = &jb_emlrtRSI;
          b_plus(&st, B_k_minus_est, r2);
        }

        emlrtMEXProfilingStatement(57, isMexOutdated);
        if (nx_tmp != i8) {
          emlrtIntegerCheckR2012b(d1, &w_emlrtDCI, (emlrtConstCTX)sp);
        }

        i8 = r7->size[0] * r7->size[1] * r7->size[2];
        r7->size[0] = 7;
        r7->size[1] = 1;
        r7->size[2] = (int8_T)(N - 1.0);
        emxEnsureCapacity_real_T(sp, r7, i8, &lc_emlrtRTEI);
        r16 = r7->data;
        for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
          r12 = _mm_loadu_pd(&S_kdot_ki_data[7 * i8 + 7 * S_kdot_ki->size[2] * j]);
          r14 = _mm_set1_pd(d2);
          _mm_storeu_pd(&r16[7 * i8], _mm_mul_pd(r14, r12));
          r12 = _mm_loadu_pd(&S_kdot_ki_data[(7 * i8 + 7 * S_kdot_ki->size[2] *
            j) + 2]);
          _mm_storeu_pd(&r16[7 * i8 + 2], _mm_mul_pd(r14, r12));
          r12 = _mm_loadu_pd(&S_kdot_ki_data[(7 * i8 + 7 * S_kdot_ki->size[2] *
            j) + 4]);
          _mm_storeu_pd(&r16[7 * i8 + 4], _mm_mul_pd(r14, r12));
          r16[7 * i8 + 6] = d2 * S_kdot_ki_data[(7 * i8 + 7 * S_kdot_ki->size[2]
            * j) + 6];
        }

        if ((r7->size[2] != S_k_est->size[2]) && ((r7->size[2] != 1) &&
             (S_k_est->size[2] != 1))) {
          emlrtDimSizeImpxCheckR2021b(r7->size[2], S_k_est->size[2], &q_emlrtECI,
            (emlrtConstCTX)sp);
        }

        if (r7->size[2] == S_k_est->size[2]) {
          d_loop_ub = 7 * r7->size[2];
          i8 = S_k_est->size[0] * S_k_est->size[1] * S_k_est->size[2];
          S_k_est->size[0] = 7;
          S_k_est->size[1] = 1;
          S_k_est->size[2] = r7->size[2];
          emxEnsureCapacity_real_T(sp, S_k_est, i8, &mc_emlrtRTEI);
          x_est_data = S_k_est->data;
          n = (d_loop_ub / 2) << 1;
          b_vectorUB = n - 2;
          for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
            r12 = _mm_loadu_pd(&r16[i8]);
            r14 = _mm_loadu_pd(&x_est_data[i8]);
            _mm_storeu_pd(&x_est_data[i8], _mm_add_pd(r12, r14));
          }

          for (i8 = n; i8 < d_loop_ub; i8++) {
            x_est_data[i8] += r16[i8];
          }
        } else {
          st.site = &ib_emlrtRSI;
          plus(&st, S_k_est, r7);
        }

        emlrtMEXProfilingStatement(58, isMexOutdated);
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b((emlrtConstCTX)sp);
        }
      }

      emlrtMEXProfilingStatement(59, isMexOutdated);
      if (r5->size[1] == r6->size[1]) {
        iv1[2] = r8->size[1];
        nx = r9->size[1];
        i8 = b_t_k->size[0] * b_t_k->size[1];
        b_t_k->size[0] = 1;
        b_t_k->size[1] = t_k->size[1];
        emxEnsureCapacity_real_T(sp, b_t_k, i8, &dc_emlrtRTEI);
        x_est_data = b_t_k->data;
        for (i8 = 0; i8 <= b_vectorUB_tmp; i8 += 2) {
          r12 = _mm_loadu_pd(&t_k_data[i8]);
          _mm_storeu_pd(&x_est_data[i8], _mm_add_pd(r12, _mm_set1_pd(b_dt_sub)));
        }

        for (i8 = b_scalarLB_tmp; i8 < i_loop_ub_tmp; i8++) {
          x_est_data[i8] = t_k_data[i8] + b_dt_sub;
        }

        i8 = r11->size[0] * r11->size[1];
        r11->size[0] = 2;
        r11->size[1] = r5->size[1];
        emxEnsureCapacity_real_T(sp, r11, i8, &rb_emlrtRTEI);
        r16 = r11->data;
        i_loop_ub_tmp = r5->size[1] << 1;
        n = (i_loop_ub_tmp / 2) << 1;
        b_vectorUB = n - 2;
        for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
          r12 = _mm_loadu_pd(&r17[i8]);
          r14 = _mm_loadu_pd(&r13[i8]);
          _mm_storeu_pd(&r16[i8], _mm_add_pd(r12, r14));
        }

        for (i8 = n; i8 < i_loop_ub_tmp; i8++) {
          r16[i8] = r17[i8] + r13[i8];
        }

        r23 = *r8;
        iv4[0] = 1;
        iv4[1] = 1;
        iv4[2] = iv1[2];
        r23.size = &iv4[0];
        r23.numDimensions = 3;
        r24 = *r9;
        iv5[0] = 1;
        iv5[1] = 1;
        iv5[2] = nx;
        r24.size = &iv5[0];
        r24.numDimensions = 3;
        st.site = &g_emlrtRSI;
        STM_diff_eq_FOH(&st, b_t_k, x_est, r11, &r23, &r24, A_k_est,
                        B_k_plus_est, B_k_minus_est, S_k_est, L, b_I, alpha, r,
                        r1, r2, r3, r4);
        r13 = r4->data;
        r16 = r3->data;
        r19 = r2->data;
        r21 = r1->data;
        r22 = r->data;
      } else {
        st.site = &g_emlrtRSI;
        binary_expand_op(&st, g_emlrtRSI, t_k, dv, b_i, dt_sub, x_est, r5, r6,
                         r8, r9, A_k_est, B_k_plus_est, B_k_minus_est, S_k_est,
                         L, b_I, alpha, r, r1, r2, r3, r4);
        r13 = r4->data;
        r16 = r3->data;
        r19 = r2->data;
        r21 = r1->data;
        r22 = r->data;
      }

      iv[0] = 7;
      iv[1] = (int32_T)(N - 1.0);
      emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &r->size[0], 2, &p_emlrtECI,
        (emlrtCTX)sp);
      for (i8 = 0; i8 < loop_ub_tmp; i8++) {
        for (nx_tmp = 0; nx_tmp < 7; nx_tmp++) {
          i9 = nx_tmp + 7 * i8;
          xdot_ki_data[i9 + 7 * xdot_ki->size[1] * (b_i + 1)] = r22[i9];
        }
      }

      iv1[0] = 7;
      iv1[1] = 7;
      iv1[2] = (int8_T)(N - 1.0);
      emlrtSubAssignSizeCheckR2012b(&iv1[0], 3, &r1->size[0], 3, &o_emlrtECI,
        (emlrtCTX)sp);
      for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
        for (nx_tmp = 0; nx_tmp < 7; nx_tmp++) {
          for (i9 = 0; i9 < 7; i9++) {
            n = (i9 + 7 * nx_tmp) + 49 * i8;
            A_kdot_ki_data[n + 49 * A_kdot_ki->size[2] * (b_i + 1)] = r21[n];
          }
        }
      }

      iv1[0] = 7;
      iv1[1] = 2;
      iv1[2] = (int8_T)(N - 1.0);
      emlrtSubAssignSizeCheckR2012b(&iv1[0], 3, &r2->size[0], 3, &n_emlrtECI,
        (emlrtCTX)sp);
      for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
        for (nx_tmp = 0; nx_tmp < 2; nx_tmp++) {
          for (i9 = 0; i9 < 7; i9++) {
            n = (i9 + 7 * nx_tmp) + 14 * i8;
            B_k_plusdot_ki_data[n + 14 * B_k_plusdot_ki->size[2] * (b_i + 1)] =
              r19[n];
          }
        }
      }

      iv1[0] = 7;
      iv1[1] = 2;
      iv1[2] = (int8_T)(N - 1.0);
      emlrtSubAssignSizeCheckR2012b(&iv1[0], 3, &r3->size[0], 3, &m_emlrtECI,
        (emlrtCTX)sp);
      for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
        for (nx_tmp = 0; nx_tmp < 2; nx_tmp++) {
          for (i9 = 0; i9 < 7; i9++) {
            n = (i9 + 7 * nx_tmp) + 14 * i8;
            B_k_minusdot_ki_data[n + 14 * B_k_minusdot_ki->size[2] * (b_i + 1)] =
              r16[n];
          }
        }
      }

      iv1[0] = 7;
      iv1[1] = 1;
      iv1[2] = (int8_T)(N - 1.0);
      emlrtSubAssignSizeCheckR2012b(&iv1[0], 3, &r4->size[0], 3, &l_emlrtECI,
        (emlrtCTX)sp);
      for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
        for (nx_tmp = 0; nx_tmp < 7; nx_tmp++) {
          i9 = nx_tmp + 7 * i8;
          S_kdot_ki_data[i9 + 7 * S_kdot_ki->size[2] * (b_i + 1)] = r13[i9];
        }
      }

      emlrtMEXProfilingStatement(60, isMexOutdated);
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b((emlrtConstCTX)sp);
      }
    }

    emlrtMEXProfilingStatement(61, isMexOutdated);
    emlrtMEXProfilingStatement(62, isMexOutdated);
    emlrtMEXProfilingStatement(63, isMexOutdated);
    emlrtMEXProfilingStatement(64, isMexOutdated);
    emlrtMEXProfilingStatement(65, isMexOutdated);
    emlrtMEXProfilingStatement(66, isMexOutdated);
    for (j = 0; j < 9; j++) {
      emlrtMEXProfilingStatement(67, isMexOutdated);
      b_dt_sub = dt_sub * dv1[j];
      i8 = r->size[0] * r->size[1];
      r->size[0] = 7;
      r->size[1] = (int32_T)(N - 1.0);
      emxEnsureCapacity_real_T(sp, r, i8, &ib_emlrtRTEI);
      r22 = r->data;
      for (i8 = 0; i8 < loop_ub_tmp; i8++) {
        r12 = _mm_loadu_pd(&xdot_ki_data[7 * i8 + 7 * xdot_ki->size[1] * j]);
        r14 = _mm_set1_pd(b_dt_sub);
        _mm_storeu_pd(&r22[7 * i8], _mm_mul_pd(r14, r12));
        r12 = _mm_loadu_pd(&xdot_ki_data[(7 * i8 + 7 * xdot_ki->size[1] * j) + 2]);
        _mm_storeu_pd(&r22[7 * i8 + 2], _mm_mul_pd(r14, r12));
        r12 = _mm_loadu_pd(&xdot_ki_data[(7 * i8 + 7 * xdot_ki->size[1] * j) + 4]);
        _mm_storeu_pd(&r22[7 * i8 + 4], _mm_mul_pd(r14, r12));
        r22[7 * i8 + 6] = b_dt_sub * xdot_ki_data[(7 * i8 + 7 * xdot_ki->size[1]
          * j) + 6];
      }

      if ((x_k->size[1] != r->size[1]) && ((x_k->size[1] != 1) && (r->size[1] !=
            1))) {
        emlrtDimSizeImpxCheckR2021b(x_k->size[1], r->size[1], &k_emlrtECI,
          (emlrtConstCTX)sp);
      }

      if (x_k->size[1] == r->size[1]) {
        b_loop_ub = 7 * x_k->size[1];
        i8 = x_k->size[0] * x_k->size[1];
        x_k->size[0] = 7;
        emxEnsureCapacity_real_T(sp, x_k, i8, &kb_emlrtRTEI);
        x_k_data = x_k->data;
        n = (b_loop_ub / 2) << 1;
        b_vectorUB = n - 2;
        for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
          r12 = _mm_loadu_pd(&x_k_data[i8]);
          r14 = _mm_loadu_pd(&r22[i8]);
          _mm_storeu_pd(&x_k_data[i8], _mm_add_pd(r12, r14));
        }

        for (i8 = n; i8 < b_loop_ub; i8++) {
          x_k_data[i8] += r22[i8];
        }
      } else {
        st.site = &rb_emlrtRSI;
        h_plus(&st, x_k, r);
        x_k_data = x_k->data;
      }

      emlrtMEXProfilingStatement(68, isMexOutdated);
      i8 = r1->size[0] * r1->size[1] * r1->size[2];
      r1->size[0] = 7;
      r1->size[1] = 7;
      r1->size[2] = (int8_T)(N - 1.0);
      emxEnsureCapacity_real_T(sp, r1, i8, &lb_emlrtRTEI);
      r21 = r1->data;
      for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
        r12 = _mm_set1_pd(b_dt_sub);
        for (nx_tmp = 0; nx_tmp < 7; nx_tmp++) {
          i9 = 7 * nx_tmp + 49 * i8;
          r14 = _mm_loadu_pd(&A_kdot_ki_data[i9 + 49 * A_kdot_ki->size[2] * j]);
          _mm_storeu_pd(&r21[i9], _mm_mul_pd(r12, r14));
          r14 = _mm_loadu_pd(&A_kdot_ki_data[((7 * nx_tmp + 49 * i8) + 49 *
            A_kdot_ki->size[2] * j) + 2]);
          _mm_storeu_pd(&r21[i9 + 2], _mm_mul_pd(r12, r14));
          r14 = _mm_loadu_pd(&A_kdot_ki_data[((7 * nx_tmp + 49 * i8) + 49 *
            A_kdot_ki->size[2] * j) + 4]);
          _mm_storeu_pd(&r21[i9 + 4], _mm_mul_pd(r12, r14));
          r21[i9 + 6] = b_dt_sub * A_kdot_ki_data[((7 * nx_tmp + 49 * i8) + 49 *
            A_kdot_ki->size[2] * j) + 6];
        }
      }

      if ((A_k->size[2] != r1->size[2]) && ((A_k->size[2] != 1) && (r1->size[2]
            != 1))) {
        emlrtDimSizeImpxCheckR2021b(A_k->size[2], r1->size[2], &j_emlrtECI,
          (emlrtConstCTX)sp);
      }

      if (A_k->size[2] == r1->size[2]) {
        b_loop_ub = 49 * A_k->size[2];
        i8 = A_k->size[0] * A_k->size[1] * A_k->size[2];
        A_k->size[0] = 7;
        A_k->size[1] = 7;
        emxEnsureCapacity_real_T(sp, A_k, i8, &nb_emlrtRTEI);
        A_k_data = A_k->data;
        n = (b_loop_ub / 2) << 1;
        b_vectorUB = n - 2;
        for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
          r12 = _mm_loadu_pd(&A_k_data[i8]);
          r14 = _mm_loadu_pd(&r21[i8]);
          _mm_storeu_pd(&A_k_data[i8], _mm_add_pd(r12, r14));
        }

        for (i8 = n; i8 < b_loop_ub; i8++) {
          A_k_data[i8] += r21[i8];
        }
      } else {
        st.site = &qb_emlrtRSI;
        g_plus(&st, A_k, r1);
        A_k_data = A_k->data;
      }

      emlrtMEXProfilingStatement(69, isMexOutdated);
      i8 = r2->size[0] * r2->size[1] * r2->size[2];
      r2->size[0] = 7;
      r2->size[1] = 2;
      r2->size[2] = (int8_T)(N - 1.0);
      emxEnsureCapacity_real_T(sp, r2, i8, &ob_emlrtRTEI);
      r19 = r2->data;
      for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
        r12 = _mm_set1_pd(b_dt_sub);
        for (nx_tmp = 0; nx_tmp < 2; nx_tmp++) {
          i9 = 7 * nx_tmp + 14 * i8;
          r14 = _mm_loadu_pd(&B_k_plusdot_ki_data[i9 + 14 * B_k_plusdot_ki->
                             size[2] * j]);
          _mm_storeu_pd(&r19[i9], _mm_mul_pd(r12, r14));
          r14 = _mm_loadu_pd(&B_k_plusdot_ki_data[((7 * nx_tmp + 14 * i8) + 14 *
            B_k_plusdot_ki->size[2] * j) + 2]);
          _mm_storeu_pd(&r19[i9 + 2], _mm_mul_pd(r12, r14));
          r14 = _mm_loadu_pd(&B_k_plusdot_ki_data[((7 * nx_tmp + 14 * i8) + 14 *
            B_k_plusdot_ki->size[2] * j) + 4]);
          _mm_storeu_pd(&r19[i9 + 4], _mm_mul_pd(r12, r14));
          r19[i9 + 6] = b_dt_sub * B_k_plusdot_ki_data[((7 * nx_tmp + 14 * i8) +
            14 * B_k_plusdot_ki->size[2] * j) + 6];
        }
      }

      if ((B_k_plus->size[2] != r2->size[2]) && ((B_k_plus->size[2] != 1) &&
           (r2->size[2] != 1))) {
        emlrtDimSizeImpxCheckR2021b(B_k_plus->size[2], r2->size[2], &i_emlrtECI,
          (emlrtConstCTX)sp);
      }

      if (B_k_plus->size[2] == r2->size[2]) {
        b_loop_ub = 14 * B_k_plus->size[2];
        i8 = B_k_plus->size[0] * B_k_plus->size[1] * B_k_plus->size[2];
        B_k_plus->size[0] = 7;
        B_k_plus->size[1] = 2;
        emxEnsureCapacity_real_T(sp, B_k_plus, i8, &pb_emlrtRTEI);
        B_k_plus_data = B_k_plus->data;
        n = (b_loop_ub / 2) << 1;
        b_vectorUB = n - 2;
        for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
          r12 = _mm_loadu_pd(&B_k_plus_data[i8]);
          r14 = _mm_loadu_pd(&r19[i8]);
          _mm_storeu_pd(&B_k_plus_data[i8], _mm_add_pd(r12, r14));
        }

        for (i8 = n; i8 < b_loop_ub; i8++) {
          B_k_plus_data[i8] += r19[i8];
        }
      } else {
        st.site = &pb_emlrtRSI;
        f_plus(&st, B_k_plus, r2);
        B_k_plus_data = B_k_plus->data;
      }

      emlrtMEXProfilingStatement(70, isMexOutdated);
      i8 = r2->size[0] * r2->size[1] * r2->size[2];
      r2->size[0] = 7;
      r2->size[1] = 2;
      r2->size[2] = (int8_T)(N - 1.0);
      emxEnsureCapacity_real_T(sp, r2, i8, &qb_emlrtRTEI);
      r19 = r2->data;
      for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
        r12 = _mm_set1_pd(b_dt_sub);
        for (nx_tmp = 0; nx_tmp < 2; nx_tmp++) {
          i9 = 7 * nx_tmp + 14 * i8;
          r14 = _mm_loadu_pd(&B_k_minusdot_ki_data[i9 + 14 *
                             B_k_minusdot_ki->size[2] * j]);
          _mm_storeu_pd(&r19[i9], _mm_mul_pd(r12, r14));
          r14 = _mm_loadu_pd(&B_k_minusdot_ki_data[((7 * nx_tmp + 14 * i8) + 14 *
            B_k_minusdot_ki->size[2] * j) + 2]);
          _mm_storeu_pd(&r19[i9 + 2], _mm_mul_pd(r12, r14));
          r14 = _mm_loadu_pd(&B_k_minusdot_ki_data[((7 * nx_tmp + 14 * i8) + 14 *
            B_k_minusdot_ki->size[2] * j) + 4]);
          _mm_storeu_pd(&r19[i9 + 4], _mm_mul_pd(r12, r14));
          r19[i9 + 6] = b_dt_sub * B_k_minusdot_ki_data[((7 * nx_tmp + 14 * i8)
            + 14 * B_k_minusdot_ki->size[2] * j) + 6];
        }
      }

      if ((B_k_minus->size[2] != r2->size[2]) && ((B_k_minus->size[2] != 1) &&
           (r2->size[2] != 1))) {
        emlrtDimSizeImpxCheckR2021b(B_k_minus->size[2], r2->size[2], &h_emlrtECI,
          (emlrtConstCTX)sp);
      }

      if (B_k_minus->size[2] == r2->size[2]) {
        b_loop_ub = 14 * B_k_minus->size[2];
        i8 = B_k_minus->size[0] * B_k_minus->size[1] * B_k_minus->size[2];
        B_k_minus->size[0] = 7;
        B_k_minus->size[1] = 2;
        emxEnsureCapacity_real_T(sp, B_k_minus, i8, &tb_emlrtRTEI);
        B_k_minus_data = B_k_minus->data;
        n = (b_loop_ub / 2) << 1;
        b_vectorUB = n - 2;
        for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
          r12 = _mm_loadu_pd(&B_k_minus_data[i8]);
          r14 = _mm_loadu_pd(&r19[i8]);
          _mm_storeu_pd(&B_k_minus_data[i8], _mm_add_pd(r12, r14));
        }

        for (i8 = n; i8 < b_loop_ub; i8++) {
          B_k_minus_data[i8] += r19[i8];
        }
      } else {
        st.site = &ob_emlrtRSI;
        f_plus(&st, B_k_minus, r2);
        B_k_minus_data = B_k_minus->data;
      }

      emlrtMEXProfilingStatement(71, isMexOutdated);
      i8 = r7->size[0] * r7->size[1] * r7->size[2];
      r7->size[0] = 7;
      r7->size[1] = 1;
      r7->size[2] = (int8_T)(N - 1.0);
      emxEnsureCapacity_real_T(sp, r7, i8, &vb_emlrtRTEI);
      r16 = r7->data;
      for (i8 = 0; i8 < c_loop_ub_tmp; i8++) {
        r12 = _mm_loadu_pd(&S_kdot_ki_data[7 * i8 + 7 * S_kdot_ki->size[2] * j]);
        r14 = _mm_set1_pd(b_dt_sub);
        _mm_storeu_pd(&r16[7 * i8], _mm_mul_pd(r14, r12));
        r12 = _mm_loadu_pd(&S_kdot_ki_data[(7 * i8 + 7 * S_kdot_ki->size[2] * j)
                           + 2]);
        _mm_storeu_pd(&r16[7 * i8 + 2], _mm_mul_pd(r14, r12));
        r12 = _mm_loadu_pd(&S_kdot_ki_data[(7 * i8 + 7 * S_kdot_ki->size[2] * j)
                           + 4]);
        _mm_storeu_pd(&r16[7 * i8 + 4], _mm_mul_pd(r14, r12));
        r16[7 * i8 + 6] = b_dt_sub * S_kdot_ki_data[(7 * i8 + 7 *
          S_kdot_ki->size[2] * j) + 6];
      }

      if ((S_k->size[2] != r7->size[2]) && ((S_k->size[2] != 1) && (r7->size[2]
            != 1))) {
        emlrtDimSizeImpxCheckR2021b(S_k->size[2], r7->size[2], &g_emlrtECI,
          (emlrtConstCTX)sp);
      }

      if (S_k->size[2] == r7->size[2]) {
        b_loop_ub = 7 * S_k->size[2];
        i8 = S_k->size[0] * S_k->size[1] * S_k->size[2];
        S_k->size[0] = 7;
        S_k->size[1] = 1;
        emxEnsureCapacity_real_T(sp, S_k, i8, &xb_emlrtRTEI);
        S_k_data = S_k->data;
        n = (b_loop_ub / 2) << 1;
        b_vectorUB = n - 2;
        for (i8 = 0; i8 <= b_vectorUB; i8 += 2) {
          r12 = _mm_loadu_pd(&S_k_data[i8]);
          r14 = _mm_loadu_pd(&r16[i8]);
          _mm_storeu_pd(&S_k_data[i8], _mm_add_pd(r12, r14));
        }

        for (i8 = n; i8 < b_loop_ub; i8++) {
          S_k_data[i8] += r16[i8];
        }
      } else {
        st.site = &nb_emlrtRSI;
        e_plus(&st, S_k, r7);
        S_k_data = S_k->data;
      }

      emlrtMEXProfilingStatement(72, isMexOutdated);
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b((emlrtConstCTX)sp);
      }
    }

    emlrtMEXProfilingStatement(73, isMexOutdated);
    emlrtMEXProfilingStatement(74, isMexOutdated);
    emlrtMEXProfilingStatement(75, isMexOutdated);
    emlrtMEXProfilingStatement(76, isMexOutdated);
    emlrtMEXProfilingStatement(77, isMexOutdated);
    emlrtMEXProfilingStatement(78, isMexOutdated);
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }

  emxFree_real_T(sp, &r11);
  emxFree_real_T(sp, &b_t_k);
  emxFree_real_T(sp, &r10);
  emxFree_real_T(sp, &r9);
  emxFree_real_T(sp, &r8);
  emxFree_real_T(sp, &r7);
  emxFree_real_T(sp, &r6);
  emxFree_real_T(sp, &r5);
  emxFree_real_T(sp, &r4);
  emxFree_real_T(sp, &r3);
  emxFree_real_T(sp, &r2);
  emxFree_real_T(sp, &r1);
  emxFree_real_T(sp, &r);
  emxFree_real_T(sp, &t_k);
  emxFree_real_T(sp, &B_k_minus_est);
  emxFree_real_T(sp, &B_k_plus_est);
  emxFree_real_T(sp, &A_k_est);
  emxFree_real_T(sp, &x_est);
  emxFree_real_T(sp, &S_kdot_ki);
  emxFree_real_T(sp, &B_k_minusdot_ki);
  emxFree_real_T(sp, &B_k_plusdot_ki);
  emxFree_real_T(sp, &A_kdot_ki);
  emxFree_real_T(sp, &xdot_ki);
  emxFree_real_T(sp, &t);
  emlrtMEXProfilingStatement(79, isMexOutdated);
  emlrtMEXProfilingStatement(80, isMexOutdated);
  emlrtMEXProfilingStatement(81, isMexOutdated);
  emlrtMEXProfilingStatement(82, isMexOutdated);
  st.site = &h_emlrtRSI;
  nx = 7 * x_k->size[1];
  b_st.site = &m_emlrtRSI;
  computeDimsData(&b_st, N - 1.0);
  n = 7;
  if (x_k->size[1] > 7) {
    n = x_k->size[1];
  }

  if ((int32_T)(N - 1.0) > muIntScalarMax_sint32(nx, n)) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:reshape_emptyReshapeLimit",
      "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }

  if ((int32_T)(N - 1.0) < 0) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
      "MATLAB:checkDimCommon:nonnegativeSize",
      "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }

  if (b_loop_ub_tmp != nx) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel",
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  if (N - 1.0 < 1.0) {
    i1 = 0;
  } else {
    if (N - 1.0 != i) {
      emlrtIntegerCheckR2012b(N - 1.0, &d_emlrtDCI, (emlrtConstCTX)sp);
    }

    if (((int32_T)(N - 1.0) < 1) || ((int32_T)(N - 1.0) > 15)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(N - 1.0), 1, 15, &b_emlrtBCI,
        (emlrtConstCTX)sp);
    }

    i1 = (int32_T)(N - 1.0);
  }

  st.site = &h_emlrtRSI;
  nx = 7 * i1;
  b_st.site = &m_emlrtRSI;
  computeDimsData(&b_st, N - 1.0);
  n = 7;
  if (i1 > 7) {
    n = i1;
  }

  if ((int32_T)(N - 1.0) > muIntScalarMax_sint32(nx, n)) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:reshape_emptyReshapeLimit",
      "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }

  if ((int32_T)(N - 1.0) < 0) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
      "MATLAB:checkDimCommon:nonnegativeSize",
      "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }

  if (b_loop_ub_tmp != nx) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel",
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  if (N - 1.0 < 1.0) {
    i = 0;
  } else {
    if (N - 1.0 != i) {
      emlrtIntegerCheckR2012b(N - 1.0, &c_emlrtDCI, (emlrtConstCTX)sp);
    }

    if (((int32_T)(N - 1.0) < 1) || ((int32_T)(N - 1.0) > 15)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(N - 1.0), 1, 15, &emlrtBCI,
        (emlrtConstCTX)sp);
    }

    i = (int32_T)(N - 1.0);
  }

  st.site = &i_emlrtRSI;
  nx = i << 1;
  b_st.site = &m_emlrtRSI;
  computeDimsData(&b_st, N - 1.0);
  n = 2;
  if (i > 2) {
    n = i;
  }

  if ((int32_T)(N - 1.0) > muIntScalarMax_sint32(nx, n)) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:reshape_emptyReshapeLimit",
      "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }

  if ((int32_T)(N - 1.0) < 0) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
      "MATLAB:checkDimCommon:nonnegativeSize",
      "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }

  i = (int32_T)(N - 1.0) << 1;
  if (i != nx) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel",
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  i1 = S_k_est->size[0] * S_k_est->size[1] * S_k_est->size[2];
  S_k_est->size[0] = 7;
  S_k_est->size[1] = 1;
  S_k_est->size[2] = (int32_T)(N - 1.0);
  emxEnsureCapacity_real_T(sp, S_k_est, i1, &ab_emlrtRTEI);
  x_est_data = S_k_est->data;
  for (i1 = 0; i1 < loop_ub_tmp; i1++) {
    for (i2 = 0; i2 < 7; i2++) {
      i3 = i2 + 7 * i1;
      x_est_data[i3] = x_ref[i3];
    }
  }

  emxInit_real_T(sp, &r18, 3, &fb_emlrtRTEI);
  st.site = &h_emlrtRSI;
  c_pagemtimes(&st, A_k, S_k_est, r18);
  r13 = r18->data;
  iv1[0] = 2;
  iv1[1] = 1;
  iv1[2] = (int32_T)(N - 1.0);
  for (i1 = 0; i1 < loop_ub_tmp; i1++) {
    nx = i1 << 1;
    u_ref_data[2 * i1] = u_ref[nx];
    u_ref_data[2 * i1 + 1] = u_ref[nx + 1];
  }

  emxInit_real_T(sp, &r20, 3, &fb_emlrtRTEI);
  st.site = &i_emlrtRSI;
  d_pagemtimes(&st, B_k_minus, u_ref_data, iv1, r20);
  r15 = r20->data;
  if ((r18->size[0] != r20->size[0]) && ((r18->size[0] != 1) && (r20->size[0] !=
        1))) {
    emlrtDimSizeImpxCheckR2021b(r18->size[0], r20->size[0], &f_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if ((r18->size[1] != r20->size[1]) && ((r18->size[1] != 1) && (r20->size[1] !=
        1))) {
    emlrtDimSizeImpxCheckR2021b(r18->size[1], r20->size[1], &e_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if ((r18->size[2] != r20->size[2]) && ((r18->size[2] != 1) && (r20->size[2] !=
        1))) {
    emlrtDimSizeImpxCheckR2021b(r18->size[2], r20->size[2], &d_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if (N < 2.0) {
    i1 = 0;
    i2 = 0;
  } else {
    i1 = 1;
    if (N != (int32_T)muDoubleScalarFloor(N)) {
      emlrtIntegerCheckR2012b(N, &b_emlrtDCI, (emlrtConstCTX)sp);
    }

    i2 = (int32_T)N;
  }

  st.site = &j_emlrtRSI;
  nx_tmp = i2 - i1;
  nx = nx_tmp << 1;
  b_st.site = &m_emlrtRSI;
  computeDimsData(&b_st, N - 1.0);
  n = 2;
  if (nx_tmp > 2) {
    n = i2 - i1;
  }

  if ((int32_T)(N - 1.0) > muIntScalarMax_sint32(nx, n)) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:reshape_emptyReshapeLimit",
      "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }

  if ((int32_T)(N - 1.0) < 0) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
      "MATLAB:checkDimCommon:nonnegativeSize",
      "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }

  if (i != nx) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel",
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  if ((r18->size[0] == r20->size[0]) && (r18->size[1] == r20->size[1]) &&
      (r18->size[2] == r20->size[2])) {
    b_loop_ub_tmp = r18->size[0] * r18->size[1] * r18->size[2];
    scalarLB = (b_loop_ub_tmp / 2) << 1;
    vectorUB = scalarLB - 2;
    for (i = 0; i <= vectorUB; i += 2) {
      r12 = _mm_loadu_pd(&r13[i]);
      r14 = _mm_loadu_pd(&r15[i]);
      _mm_storeu_pd(&r13[i], _mm_add_pd(r12, r14));
    }

    for (i = scalarLB; i < b_loop_ub_tmp; i++) {
      r13[i] += r15[i];
    }
  } else {
    st.site = &h_emlrtRSI;
    j_plus(&st, r18, r20);
    r13 = r18->data;
  }

  iv1[0] = 2;
  iv1[1] = 1;
  iv1[2] = (int32_T)(N - 1.0);
  for (i = 0; i < loop_ub_tmp; i++) {
    nx = (i1 + i) << 1;
    u_ref_data[2 * i] = u_ref[nx];
    u_ref_data[2 * i + 1] = u_ref[nx + 1];
  }

  st.site = &j_emlrtRSI;
  d_pagemtimes(&st, B_k_plus, u_ref_data, iv1, r20);
  r15 = r20->data;
  if ((r18->size[0] != r20->size[0]) && ((r18->size[0] != 1) && (r20->size[0] !=
        1))) {
    emlrtDimSizeImpxCheckR2021b(r18->size[0], r20->size[0], &f_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if ((r18->size[1] != r20->size[1]) && ((r18->size[1] != 1) && (r20->size[1] !=
        1))) {
    emlrtDimSizeImpxCheckR2021b(r18->size[1], r20->size[1], &e_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if ((r18->size[2] != r20->size[2]) && ((r18->size[2] != 1) && (r20->size[2] !=
        1))) {
    emlrtDimSizeImpxCheckR2021b(r18->size[2], r20->size[2], &d_emlrtECI,
      (emlrtConstCTX)sp);
  }

  st.site = &k_emlrtRSI;
  emlrtMEXProfilingFunctionEntry((char_T *)zero_if_empty_complete_name,
    isMexOutdated);

  /* ZERO_IF_EMPTY Summary of this function goes here */
  /*    Detailed explanation goes here */
  emlrtMEXProfilingStatement(2, isMexOutdated);
  emlrtMEXProfilingStatement(4, isMexOutdated);
  emlrtMEXProfilingFunctionExit(isMexOutdated);
  if ((r18->size[0] == r20->size[0]) && (r18->size[1] == r20->size[1]) &&
      (r18->size[2] == r20->size[2])) {
    b_loop_ub_tmp = r18->size[0] * r18->size[1] * r18->size[2];
    scalarLB = (b_loop_ub_tmp / 2) << 1;
    vectorUB = scalarLB - 2;
    for (i = 0; i <= vectorUB; i += 2) {
      r12 = _mm_loadu_pd(&r13[i]);
      r14 = _mm_loadu_pd(&r15[i]);
      _mm_storeu_pd(&r13[i], _mm_add_pd(r12, r14));
    }

    for (i = scalarLB; i < b_loop_ub_tmp; i++) {
      r13[i] += r15[i];
    }
  } else {
    st.site = &h_emlrtRSI;
    j_plus(&st, r18, r20);
    r13 = r18->data;
  }

  emxFree_real_T(sp, &r20);
  i = S_k_est->size[0] * S_k_est->size[1] * S_k_est->size[2];
  S_k_est->size[0] = 7;
  S_k_est->size[1] = 1;
  S_k_est->size[2] = S_k->size[2];
  emxEnsureCapacity_real_T(sp, S_k_est, i, &eb_emlrtRTEI);
  x_est_data = S_k_est->data;
  b_loop_ub_tmp = 7 * S_k->size[2];
  scalarLB = (b_loop_ub_tmp / 2) << 1;
  vectorUB = scalarLB - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r12 = _mm_loadu_pd(&S_k_data[i]);
    _mm_storeu_pd(&x_est_data[i], _mm_mul_pd(r12, _mm_set1_pd(0.0)));
  }

  for (i = scalarLB; i < b_loop_ub_tmp; i++) {
    x_est_data[i] = S_k_data[i] * 0.0;
  }

  emxInit_real_T(sp, &r25, 3, &fb_emlrtRTEI);
  st.site = &k_emlrtRSI;
  zero_if_empty(&st, S_k_est, r25);
  emxFree_real_T(sp, &S_k_est);
  if ((r18->size[0] != r25->size[0]) && ((r18->size[0] != 1) && (r25->size[0] !=
        1))) {
    emlrtDimSizeImpxCheckR2021b(r18->size[0], r25->size[0], &f_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if ((r18->size[1] != r25->size[1]) && ((r18->size[1] != 1) && (r25->size[1] !=
        1))) {
    emlrtDimSizeImpxCheckR2021b(r18->size[1], r25->size[1], &e_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if ((r18->size[2] != r25->size[2]) && ((r18->size[2] != 1) && (r25->size[2] !=
        1))) {
    emlrtDimSizeImpxCheckR2021b(r18->size[2], r25->size[2], &d_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if ((r18->size[0] == r25->size[0]) && (r18->size[1] == r25->size[1]) &&
      (r18->size[2] == r25->size[2])) {
    b_loop_ub_tmp = r18->size[0] * r18->size[1] * r18->size[2];
    i = r25->size[0] * r25->size[1] * r25->size[2];
    r25->size[0] = r18->size[0];
    r25->size[1] = r18->size[1];
    r25->size[2] = r18->size[2];
    emxEnsureCapacity_real_T(sp, r25, i, &fb_emlrtRTEI);
    r15 = r25->data;
    scalarLB = (b_loop_ub_tmp / 2) << 1;
    vectorUB = scalarLB - 2;
    for (i = 0; i <= vectorUB; i += 2) {
      r12 = _mm_loadu_pd(&r13[i]);
      r14 = _mm_loadu_pd(&r15[i]);
      _mm_storeu_pd(&r15[i], _mm_add_pd(r12, r14));
    }

    for (i = scalarLB; i < b_loop_ub_tmp; i++) {
      r15[i] += r13[i];
    }
  } else {
    st.site = &h_emlrtRSI;
    i_plus(&st, r25, r18);
    r15 = r25->data;
  }

  emxFree_real_T(sp, &r18);
  if ((r25->size[0] != 7) && (r25->size[0] != 1)) {
    emlrtDimSizeImpxCheckR2021b(7, r25->size[0], &c_emlrtECI, (emlrtConstCTX)sp);
  }

  if (((int32_T)(N - 1.0) != r25->size[2]) && (((int32_T)(N - 1.0) != 1) &&
       (r25->size[2] != 1))) {
    emlrtDimSizeImpxCheckR2021b((int32_T)(N - 1.0), r25->size[2], &b_emlrtECI,
      (emlrtConstCTX)sp);
  }

  if ((r25->size[0] == 7) && ((int32_T)(N - 1.0) == r25->size[2])) {
    i = d_k->size[0] * d_k->size[1] * d_k->size[2];
    d_k->size[0] = 7;
    loop_ub = r25->size[1];
    d_k->size[1] = r25->size[1];
    d_k->size[2] = (int32_T)(N - 1.0);
    emxEnsureCapacity_real_T(sp, d_k, i, &hb_emlrtRTEI);
    x_est_data = d_k->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        r12 = _mm_loadu_pd(&x_k_data[7 * i]);
        r14 = _mm_loadu_pd(&r15[r25->size[0] * i1 + r25->size[0] * r25->size[1] *
                           i]);
        _mm_storeu_pd(&x_est_data[7 * i1 + 7 * d_k->size[1] * i], _mm_sub_pd(r12,
          r14));
        r12 = _mm_loadu_pd(&x_k_data[7 * i + 2]);
        r14 = _mm_loadu_pd(&r15[(r25->size[0] * i1 + r25->size[0] * r25->size[1]
          * i) + 2]);
        _mm_storeu_pd(&x_est_data[(7 * i1 + 7 * d_k->size[1] * i) + 2],
                      _mm_sub_pd(r12, r14));
        r12 = _mm_loadu_pd(&x_k_data[7 * i + 4]);
        r14 = _mm_loadu_pd(&r15[(r25->size[0] * i1 + r25->size[0] * r25->size[1]
          * i) + 4]);
        _mm_storeu_pd(&x_est_data[(7 * i1 + 7 * d_k->size[1] * i) + 4],
                      _mm_sub_pd(r12, r14));
        x_est_data[(7 * i1 + 7 * d_k->size[1] * i) + 6] = x_k_data[7 * i + 6] -
          r15[(r25->size[0] * i1 + r25->size[0] * r25->size[1] * i) + 6];
      }
    }
  } else {
    st.site = &h_emlrtRSI;
    binary_expand_op_8(&st, d_k, x_k, N, r25);
  }

  emxFree_real_T(sp, &r25);
  emlrtMEXProfilingStatement(83, isMexOutdated);
  if (N < 2.0) {
    i = 0;
    i1 = 0;
  } else {
    i = 1;
    if (N != (int32_T)muDoubleScalarFloor(N)) {
      emlrtIntegerCheckR2012b(N, &emlrtDCI, (emlrtConstCTX)sp);
    }

    i1 = (int32_T)N;
  }

  loop_ub = x_k->size[1];
  i2 = i1 - i;
  if ((x_k->size[1] != i2) && ((x_k->size[1] != 1) && (i2 != 1))) {
    emlrtDimSizeImpxCheckR2021b(x_k->size[1], i2, &emlrtECI, (emlrtConstCTX)sp);
  }

  if (x_k->size[1] == i2) {
    i1 = x_k->size[0] * x_k->size[1];
    x_k->size[0] = 7;
    emxEnsureCapacity_real_T(sp, x_k, i1, &jb_emlrtRTEI);
    x_k_data = x_k->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      r12 = _mm_loadu_pd(&x_k_data[7 * i1]);
      i2 = 7 * (i + i1);
      _mm_storeu_pd(&x_k_data[7 * i1], _mm_sub_pd(r12, _mm_loadu_pd(&x_ref[i2])));
      i3 = 7 * i1 + 2;
      r12 = _mm_loadu_pd(&x_k_data[i3]);
      _mm_storeu_pd(&x_k_data[i3], _mm_sub_pd(r12, _mm_loadu_pd(&x_ref[i2 + 2])));
      i3 = 7 * i1 + 4;
      r12 = _mm_loadu_pd(&x_k_data[i3]);
      _mm_storeu_pd(&x_k_data[i3], _mm_sub_pd(r12, _mm_loadu_pd(&x_ref[i2 + 4])));
      i3 = 7 * i1 + 6;
      x_k_data[i3] -= x_ref[i2 + 6];
    }

    st.site = &l_emlrtRSI;
    vecnorm(&st, x_k, Delta);
  } else {
    st.site = &l_emlrtRSI;
    binary_expand_op_7(&st, Delta, l_emlrtRSI, x_k, x_ref, i, i1);
  }

  emxFree_real_T(sp, &x_k);
  emlrtMEXProfilingStatement(84, isMexOutdated);
  emlrtMEXProfilingFunctionExit(isMexOutdated);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (discretize_error_dynamics_FOH_RKV65_3DoF.c) */
