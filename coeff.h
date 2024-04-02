/* vim: set tabstop=4:softtabstop=4:shiftwidth=4:noexpandtab */

/*

2.4 kbps MELP Proposed Federal Standard speech coder

*/

/* coeff.h: filter coefficient header file */
/*        定义常数    h文件声明变量          */
/*                                         */

#ifndef _COEFF_H_
#define _COEFF_H_

#include "sc1200.h"

/* Lowpass filter coefficient in second-order sections */

extern const int16_t lpf_num[];
extern const int16_t lpf_den[];

/* Butterworth bandpass filters in second-order sections */
extern const int16_t bpf_num[];
extern const int16_t bpf_num_class[];

/* sign of coefficients for bpf_den is reversed */
extern const int16_t bpf_den[];
extern const int16_t bpf_den_class[];

/* Hamming window coefficents in Q15 */
extern const int16_t win_cof[];

/* Bandpass filter coeffients */
extern const int16_t bp_cof[][MIX_ORD + 1];

/* Triangle pulse dispersion filter */
extern const int16_t disp_cof[];

#endif
