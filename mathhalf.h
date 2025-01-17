/* vim: set tabstop=4:softtabstop=4:shiftwidth=4:noexpandtab */

/*
2.4 kbps MELP Proposed Federal Standard speech coder
*/

#ifndef _MATHHALF_H_
#define _MATHHALF_H_

/* addition */

static inline int16_t melpe_add(int16_t var1, int16_t var2) __attribute__((always_inline));	/* 1 ops */

static inline int16_t melpe_sub(int16_t var1, int16_t var2) __attribute__((always_inline));	/* 1 ops */

static inline int32_t melpe_L_add(int32_t L_var1, int32_t L_var2) __attribute__((always_inline));	/* 2 ops */

static inline int32_t melpe_L_sub(int32_t L_var1, int32_t L_var2) __attribute__((always_inline));	/* 2 ops */

/* multiplication */

static inline int16_t melpe_mult(int16_t var1, int16_t var2) __attribute__((always_inline));	/* 1 ops */

static inline int32_t melpe_L_mult(int16_t var1, int16_t var2) __attribute__((always_inline));	/* 1 ops */

/* arithmetic shifts */

static inline int16_t melpe_shr(int16_t var1, int16_t var2);	/* 1 ops */

static inline int16_t melpe_shl(int16_t var1, int16_t var2) __attribute__((always_inline));	/* 1 ops */

static inline int32_t melpe_L_shr(int32_t L_var1, int16_t var2);	/* 2 ops */

static inline int32_t melpe_L_shl(int32_t L_var1, int16_t var2) __attribute__((always_inline));	/* 2 ops */

static inline int16_t melpe_shift_r(int16_t var, int16_t var2) __attribute__((always_inline));	/* 2 ops */

static inline int32_t melpe_L_shift_r(int32_t L_var, int16_t var2) __attribute__((always_inline));	/* 3 ops */

/* absolute value  */

static inline int16_t melpe_abs_s(int16_t var1) __attribute__((always_inline));	/* 1 ops */

static inline int32_t melpe_L_abs(int32_t var1) __attribute__((always_inline));	/* 3 ops */

/* multiply accumulate	*/

static inline int32_t melpe_L_mac(int32_t L_var3, int16_t var1, int16_t var2) __attribute__((always_inline));	/* 1 op */

static inline int32_t melpe_L_msu(int32_t L_var3, int16_t var1, int16_t var2) __attribute__((always_inline));	/* 1 op */

static inline int16_t melpe_msu_r(int32_t L_var3, int16_t var1, int16_t var2) __attribute__((always_inline));	/* 2 op */

/* negation  */

static inline int16_t melpe_negate(int16_t var1) __attribute__((always_inline));	/* 1 ops */

static inline int32_t melpe_L_negate(int32_t L_var1) __attribute__((always_inline));	/* 2 ops */

/* Accumulator manipulation */

static inline int32_t melpe_L_deposit_l(int16_t var1) __attribute__((always_inline));	/* 1 ops */

static inline int32_t melpe_L_deposit_h(int16_t var1) __attribute__((always_inline));	/* 1 ops */

static inline int16_t melpe_extract_l(int32_t L_var1) __attribute__((always_inline));	/* 1 ops */

static inline int16_t melpe_extract_h(int32_t L_var1) __attribute__((always_inline));	/* 1 ops */

/* r_ound */

static inline int16_t melpe_r_ound(int32_t L_var1) __attribute__((always_inline));	/* 1 ops */

/* Normalization */

static inline int16_t melpe_norm_l(int32_t L_var1) __attribute__((always_inline));	/* 30 ops */

static inline int16_t melpe_norm_s(int16_t var1) __attribute__((always_inline));	/* 15 ops */

/* Division */

static inline int16_t melpe_divide_s(int16_t var1, int16_t var2) __attribute__((always_inline));	/* 18 ops */

/* -------------------------------------------------------------------------- */
/* 40-Bit Routines....added by Andre 11/23/99 */

/* new 40 bits basic operators */

static inline int64_t melpe_L40_add(int64_t acc, int32_t L_var1) __attribute__((always_inline));

static inline int64_t melpe_L40_sub(int64_t acc, int32_t L_var1) __attribute__((always_inline));

static inline int64_t melpe_L40_mac(int64_t acc, int16_t var1, int16_t var2) __attribute__((always_inline));

static inline int64_t melpe_L40_msu(int64_t acc, int16_t var1, int16_t var2) __attribute__((always_inline));

static inline int64_t melpe_L40_shl(int64_t acc, int16_t var1);

static inline int64_t melpe_L40_shr(int64_t acc, int16_t var1) __attribute__((always_inline));

static inline int64_t melpe_L40_negate(int64_t acc) __attribute__((always_inline));

static inline int16_t melpe_norm32(int64_t acc) __attribute__((always_inline));
static inline int32_t melpe_L_sat32(int64_t acc) __attribute__((always_inline));

#include "mathhalf_i.h"

#endif
