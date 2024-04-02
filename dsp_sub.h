/* vim: set tabstop=4:softtabstop=4:shiftwidth=4:noexpandtab */

/*	2.4 kbps MELP Proposed Federal Standard speech coder	*/
/* ======================= */
/* dsp_sub.h: include file */
/* ======================= */

#ifndef _DSP_SUB_H_
#define _DSP_SUB_H_

#include <stdio.h>

void envelope(int16_t input[], int16_t prev_in, int16_t output[],	//	计算输入信号的包络
	      int16_t npts);

void fill(int16_t output[], int16_t fillval, int16_t npts);		// 填充给定value

void L_fill(int32_t output[], int32_t fillval, int16_t npts);		// 填充的是 long

void interp_array(int16_t prev[], int16_t curr[], int16_t out[],	// 对数 插值
		  int16_t ifact, int16_t size);

int16_t median3(int16_t input[]);	// 计算输入数组中的元素中值，与其前后两个元素相关

void pack_code(int16_t code, unsigned char **ptr_ch_begin,
	       int16_t * ptr_ch_bit, int16_t numbits, int16_t wsize);	//	将一个整数编码为一系列位，将编码结果存储在字节数组中

int16_t peakiness(int16_t input[], int16_t npts);	// 计算输入数组的峰值指标

void quant_u(int16_t * p_data, int16_t * p_index, int16_t qmin,	// 对输入数据进行量化
	     int16_t qmax, int16_t nlev, int16_t nlev_q,
	     int16_t double_flag, int16_t scale);

void quant_u_dec(int16_t index, int16_t * p_data, int16_t qmin,	// 对量化后的数据进行解码，还原为原始数据
		 int16_t qmax, int16_t nlev_q, int16_t scale);

void rand_num(int16_t output[], int16_t amplitude, int16_t npts);	//	生成具有指定幅度随机数数组

int16_t rand_minstdgen();	// 生成基于伪随机算法minstdgen的伪随机数

BOOLEAN unpack_code(unsigned char **ptr_ch_begin, int16_t * ptr_ch_bit,	//	从字节数组中解码指定位数的编码，解码结果存储在整数变量中
		    int16_t * code, int16_t numbits, int16_t wsize,
		    uint16_t erase_mask);

void window(int16_t input[], const int16_t win_coeff[],	// 对输入信号应用窗函数
	    int16_t output[], int16_t npts);

void window_Q(int16_t input[], int16_t win_coeff[], int16_t output[],	// 在上述函数上添加指定的Q指值
	      int16_t npts, int16_t Qin);

void zerflt(int16_t input[], const int16_t coeff[], int16_t output[],	// 0 相位滤波器，对输入信号滤波
	    int16_t order, int16_t npts);

void zerflt_Q(int16_t input[], const int16_t coeff[],	// 0 相位滤波器，但 Q 值被指定
	      int16_t output[], int16_t order, int16_t npts,
	      int16_t Q_coeff);

void iir_2nd_d(int16_t input[], const int16_t den[],	// 二阶直接IIR滤波器对输入信号滤波
	       const int16_t num[], int16_t output[], int16_t delin[],
	       int16_t delout_hi[], int16_t delout_lo[], int16_t npts);

void iir_2nd_s(int16_t input[], const int16_t den[],	// 二阶级联IIR滤波器对输入信号滤波
	       const int16_t num[], int16_t output[],
	       int16_t delin[], int16_t delout[], int16_t npts);

int16_t interp_scalar(int16_t prev, int16_t curr, int16_t ifact);	// 对两个数进行插值，根据插值因子，两个数平均值插值进去
//	返回的是一个整数
#endif
