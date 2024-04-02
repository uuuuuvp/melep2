/* vim: set tabstop=4:softtabstop=4:shiftwidth=4:noexpandtab */

/* ================================================================== */
/*                                                                    */
/*    Microsoft Speech coder     ANSI-C Source Code                   */
/*    SC1200 1200 bps speech coder                                    */
/*    Fixed Point Implementation      Version 7.0                     */
/*    Copyright (C) 2000, Microsoft Corp.                             */
/*    All rights reserved.                                            */
/*                                                                    */
/* ================================================================== */

/*		2.4 kbps MELP Proposed Federal Standard speech coder		*/

#include "sc1200.h"
#include "mathhalf.h"
#include "macro.h"
#include "lpc_lib.h"
#include "mat_lib.h"
#include "vq_lib.h"
#include "fs_lib.h"
#include "math_lib.h"
#include "constant.h"
#include "global.h"
#include "harm.h"
#include "fsvq_cb.h"
#include "dsp_sub.h"
#include "melp_sub.h"
#include "coeff.h"

#if POSTFILTER
#include "postfilt.h"
#endif

#define INV_LPC_ORD			2979
			       /* ((1.0/(LPC_ORD + 1))*(1<<15) + 0.5) for Q15 */
#define X005_Q19			26214	/* 0.05 * (1 << 19) */
#define X025_Q15			8192	/* 0.25 * (1 << 15) */
#define X1_732_Q14			28377	/* 1.732 * (1 << 14) */
#define X12_Q8				3072	/* 12.0 * (1 << 8) */
#define X30_Q8				7680	/* 30.0 * (1 << 8) */
#define TILT_ORD			1
#define SCALEOVER			10
#define INV_SCALEOVER_Q18	26214	/* ((1.0/SCALEOVER)*(1 << 18)) */
#define PDEL				SCALEOVER

#if (MIX_ORD > DISP_ORD)
#define BEGIN MIX_ORD
#else
#define BEGIN DISP_ORD		// 64
#endif

#define ORIGINAL_SYNTH_GAIN		FALSE
#if ORIGINAL_SYNTH_GAIN
#define SYN_GAIN_Q4				16000	/* (1000.0*(1 << 4)) */
#else
#define SYN_GAIN_Q4				32000
#endif

static struct melp_param prev_par;
static int16_t sigsave[PITCHMAX];
static int16_t syn_begin;
static BOOLEAN erase;

/* Prototype */

static void melp_syn(struct melp_param *par, int16_t sp_out[]);

void rmquan(struct melp_param *par, int16_t sp_out[], struct melp_param *par2)		// 在 sc1200.h 里面定义函数
{
	register int16_t i;
		/* Copy previous period of processed speech to output array */
	if (syn_begin > 0) {
		if (syn_begin > frameSize) {	/* impossible */
			v_equ(sp_out, sigsave, frameSize);
			/* past end: save remainder in sigsave[0] */
			v_equ(sigsave, &sigsave[frameSize],
			      (int16_t) (syn_begin - frameSize));	// framesize 在 melp_syn 是在 melp_init 里，1200 初始化为 540
		} else
			v_equ(sp_out, sigsave, syn_begin);
	}
	erase = FALSE;		/* no erasures yet */	// 表示读取到的数据正常，不用擦除，继承前一帧的数据
// 相当于是给 erase 做初始化
	/* Read and decode channel input buffer. */
// 	if (rate == RATE2400)
// 		erase = melp_chn_read(&quant_par, par, &prev_par, chbuf);
// #if !SKIP_CHANNEL
// 	else
// 		erase = (BOOLEAN) low_rate_chn_read(&quant_par, par, &prev_par);
// #endif
	*par = *par2;
	if (rate == RATE2400) {
		// par->uv_flag = quant_par.uv_flag[0];
		par->uv_flag = par2->uv_flag;
		melp_syn(par2, sp_out);
	} else {
		for (i = 0; i < NF; i++) {
			melp_syn(&par[i], &sp_out[i * FRAME]);
			if ((syn_begin > 0) && (i < NF - 1))
				v_equ(&sp_out[(i + 1) * FRAME], sigsave,
				      syn_begin);
		}
	}
}

/****************************************************************************
**
** Function:		synthesis
**
** Description: 	The synthesis routine for the sc1200 coder
**
** Arguments:
**
**	melp_param *par ---- output encoded melp parameters
**	int16_t sp_in[] ---- input speech data buffer
**
** Return value:	None
**
*****************************************************************************/
void synthesis(struct melp_param *par, int16_t sp_out[])
{
	register int16_t i;

	/* Copy previous period of processed speech to output array */
	if (syn_begin > 0) {
		if (syn_begin > frameSize) {	/* impossible */
			v_equ(sp_out, sigsave, frameSize);
			/* past end: save remainder in sigsave[0] */
			v_equ(sigsave, &sigsave[frameSize],
			      (int16_t) (syn_begin - frameSize));	// framesize 在 melp_syn 是在 melp_init 里，1200 初始化为 540
		} else
			v_equ(sp_out, sigsave, syn_begin);
	}

	erase = FALSE;		/* no erasures yet */	// 表示读取到的数据正常，不用擦除，继承前一帧的数据
// 相当于是给 erase 做初始化
	/* Read and decode channel input buffer. */
	if (rate == RATE2400)
		erase = melp_chn_read(&quant_par, par, &prev_par, chbuf);
#if !SKIP_CHANNEL
	else
		erase = (BOOLEAN) low_rate_chn_read(&quant_par, par, &prev_par);
#endif

	if (rate == RATE2400) {
		// par->uv_flag = quant_par.uv_flag[0];
		melp_syn(par, sp_out);
	} else {
		for (i = 0; i < NF; i++) {
			melp_syn(&par[i], &sp_out[i * FRAME]);
			if ((syn_begin > 0) && (i < NF - 1))
				v_equ(&sp_out[(i + 1) * FRAME], sigsave,
				      syn_begin);
		}
	}
}

/* Name: melp_syn.c                                                           */
/*  Description: MELP synthesis                                               */
/*    This program takes the new parameters for a speech                      */
/*    frame and synthesizes the output speech.  It keeps                      */
/*    an internal record of the previous frame parameters                     */
/*    to use for interpolation.                                               */
/*  Inputs:                                                                   */
/*    *par - MELP parameter structure                                         */
/*  Outputs:                                                                  */
/*    speech[] - output speech signal                                         */
/*  Returns: void                                                             */

static void melp_syn(struct melp_param *par, int16_t sp_out[])
{
	register int16_t i;
	static BOOLEAN firstTime = TRUE;
	static int16_t noise_gain = MIN_NOISE_Q8;
	static int16_t prev_lpc_gain = ONE_Q15;
	static int16_t lpc_del[LPC_ORD];	/* Q0 */
	static int16_t prev_tilt;	// 频谱倾斜程度
	static int16_t prev_pcof[MIX_ORD + 1], prev_ncof[MIX_ORD + 1];	// 周期以及噪声系数
	static int16_t disp_del[DISP_ORD];		// 64 大小
	static int16_t ase_del[LPC_ORD], tilt_del[TILT_ORD];
	static int16_t pulse_del[MIX_ORD], noise_del[MIX_ORD];
	int16_t fs_real[PITCHMAX];		// PITCHMAX 160
	int16_t sig2[BEGIN + PITCHMAX];
	int16_t sigbuf[BEGIN + PITCHMAX];
	int16_t gaincnt, length;
	int16_t intfact, intfact1, ifact, ifact_gain;	// ifact 为插值因子
	int16_t gain, pulse_gain, pitch, jitter;
	int16_t curr_tilt, tilt_cof[TILT_ORD + 1];
	int16_t sig_prob, syn_gain, lpc_gain;
	int16_t lsf[LPC_ORD];
	int16_t lpc[LPC_ORD + 1];
	int16_t ase_num[LPC_ORD + 1], ase_den[LPC_ORD];
	int16_t curr_pcof[MIX_ORD + 1], curr_ncof[MIX_ORD + 1];
	int16_t pulse_cof[MIX_ORD + 1], noise_cof[MIX_ORD + 1];
	int16_t temp1, temp2;
	int32_t L_temp1, L_temp2;
	int16_t fc_prev, fc_curr, fc;

	/* Update adaptive noise level estimate based on last gain */
	if (firstTime) {	// 当 first 时候，2400 每帧都是第一次， 1200 只有首帧
		noise_gain = par->gain[NUM_GAINFR - 1];	/* noise_gain in Q8 */
		prev_tilt = 0;
		v_zap(prev_pcof, MIX_ORD + 1);
		v_zap(prev_ncof, MIX_ORD + 1);
		prev_ncof[MIX_ORD / 2] = ONE_Q15;
		v_zap(disp_del, DISP_ORD);
		v_zap(ase_del, LPC_ORD);
		v_zap(tilt_del, TILT_ORD);
		v_zap(pulse_del, MIX_ORD);
		v_zap(noise_del, MIX_ORD);
		firstTime = FALSE;	// first 初始化已完成
	} else if (!erase) {	// 当 erase 为真，表示当前帧语音信号被擦除
		for (i = 0; i < NUM_GAINFR; i++) {	// 循环增益 ，增益数组就这么大
			noise_est(par->gain[i], &noise_gain, UPCONST_Q19,
				  DOWNCONST_Q17, MIN_NOISE_Q8, MAX_NOISE_Q8);

			/* Adjust gain based on noise level (noise suppression) */
			noise_sup(&par->gain[i], noise_gain, MAX_NS_SUP_Q8,
				  MAX_NS_ATT_Q8, NFACT_Q8);
		}
	}	// if (!erase) {专为 1200；因为2400 每帧都得调用 melp_syn

	if (par->uv_flag && (rate == RATE1200)) {	// 清浊音帧判断
		fill(par->fs_mag, ONE_Q13, NUM_HARM);
		par->pitch = UV_PITCH_Q7;	// 6400
		par->jitter = X025_Q15;		// 0.25*2^15
	}
// 以上多为 初始化， 赋个值等等
	/* Un-weight Fourier magnitudes */
	if (!par->uv_flag && !erase)
		window_Q(par->fs_mag, w_fs_inv, par->fs_mag, NUM_HARM, 14);
		// 加窗，，编码那一步也有加窗
	/* Clamp LSP bandwidths to avoid sharp LPC filters */
	lpc_clamp(par->lsf, BWMIN_Q15, LPC_ORD);

	/* Calculate spectral tilt for current frame for spectral enhancement */
	// 计算频谱倾斜程度， 用于频谱增强
	tilt_cof[0] = ONE_Q15;	/* tilt_cof in Q15 */
	// 将 lsp 带宽 限制到某个最小值
	lpc_lsp2pred(par->lsf, &(lpc[1]), LPC_ORD);

	/* Use LPC prediction gain for adaptive scaling */
	// 使用线性预测模型进行自适应缩放
	/*      lpc_gain = sqrt(lpc_pred2refl(lpc, sig2, LPC_ORD)); */
	/* Here we only make use of sig2[0] from the returned value instead of    */
	/* using the whole array sig2[].                                          */
	// 计算 lpc 增益
	lpc_gain = lpc_pred2refl(&(lpc[1]), sig2, LPC_ORD);
	// 计算 增益平方根
	lpc_gain = sqrt_fxp(lpc_gain, 15);	/* lpc_gain in Q15 */
	// 频谱倾斜程度
	if (sig2[0] < 0)
		curr_tilt = melpe_shr(sig2[0], 1);	/* curr_tilt in Q15 */
	else
		curr_tilt = 0;

	/* Disable pitch interpolation for high-pitched onsets */

	/*      if (par->pitch < 0.5*prev_par.pitch &&
	   par->gain[0] > 6.0 + prev_par.gain[NUM_GAINFR - 1]) */

	temp1 = melpe_shr(prev_par.pitch, 1);
	temp2 = melpe_add(SIX_Q8, prev_par.gain[NUM_GAINFR - 1]);
	if ((par->pitch < temp1) && (par->gain[0] > temp2)) {
		/* copy current pitch into previous */
		prev_par.pitch = par->pitch;
	}

	/* Set pulse and noise coefficients based on voicing strengths */
	v_zap(curr_pcof, MIX_ORD + 1);	/* curr_pcof in Q14 */	
	v_zap(curr_ncof, MIX_ORD + 1);	/* curr_ncof in Q14 */
	for (i = 0; i < NUM_BANDS; i++) {
		if (par->bpvc[i] > X05_Q14)
			v_add(curr_pcof, bp_cof[i], MIX_ORD + 1);	/* bp_cof in Q14 */
		else
			v_add(curr_ncof, bp_cof[i], MIX_ORD + 1);
	}

	/* Process each pitch period */
	while (syn_begin < FRAME) {

		/* interpolate previous and current parameters */
		ifact = melpe_divide_s(syn_begin, FRAME);	/* ifact in Q15 */
		// 计算插值因子
		if (syn_begin >= GAINFR) {
			gaincnt = 2;
			temp1 = melpe_sub(syn_begin, GAINFR);
			ifact_gain = melpe_divide_s(temp1, GAINFR);	// GAINFR 为常数
		} else {
			gaincnt = 1;
			ifact_gain = melpe_divide_s(syn_begin, GAINFR);
		}

		/* interpolate gain.  It is assumed that par->gain[] are obtained     */
		/* from gain_vq_cb[] in "qnt12_cb.c", and gain_vq_cb[] lies between   */
		/* 2564 and 18965 (Q8).  Therefore, the interpolated "gain" is also   */
		/* assumed to be between these two values.                            */
		if (gaincnt > 1) {
			/*      gain = ifact_gain * par->gain[gaincnt - 1] +
			   (1.0 - ifact_gain) * par->gain[gaincnt - 2];               */
			L_temp1 = melpe_L_mult(par->gain[gaincnt - 1], ifact_gain);
			temp1 = melpe_sub(ONE_Q15, ifact_gain);
			L_temp2 = melpe_L_mult(par->gain[gaincnt - 2], temp1);
			gain = melpe_extract_h(melpe_L_add(L_temp1, L_temp2));	/* gain in Q8 */
		} else {
			/*      gain = ifact_gain * par->gain[gaincnt - 1] +
			   (1.0 - ifact_gain) * prev_par.gain[NUM_GAINFR - 1];        */
			L_temp1 = melpe_L_mult(par->gain[gaincnt - 1], ifact_gain);
			temp1 = melpe_sub(ONE_Q15, ifact_gain);
			L_temp2 = melpe_L_mult(prev_par.gain[NUM_GAINFR - 1], temp1);
			gain = melpe_extract_h(melpe_L_add(L_temp1, L_temp2));	/* gain in Q8 */
		}

/* Set overall interpolation path based on gain change */

		temp1 =
		    melpe_sub(par->gain[NUM_GAINFR - 1],
			prev_par.gain[NUM_GAINFR - 1]);
		if (melpe_abs_s(temp1) > SIX_Q8) {
			/* Power surge: use gain adjusted interpolation */
			/*      intfact = (gain - prev_par.gain[NUM_GAINFR - 1])/temp; */
			temp2 = melpe_sub(gain, prev_par.gain[NUM_GAINFR - 1]);
			if (((temp2 > 0) && (temp1 < 0)) ||
			    ((temp2 < 0) && (temp1 > 0)))
				intfact = 0;
			else {
				temp1 = melpe_abs_s(temp1);
				temp2 = melpe_abs_s(temp2);
				if (temp2 >= temp1)
					intfact = ONE_Q15;
				else
					intfact = melpe_divide_s(temp2, temp1);	/* intfact in Q15 */
			}
		} else		/* Otherwise, linear interpolation */
			intfact = ifact;

		/* interpolate LSF's and convert to LPC filter */
		interp_array(prev_par.lsf, par->lsf, lsf, intfact, LPC_ORD);
		//这行代码调用interp_array函数来进行插值操作。它使用prev_par.lsf（前一个帧的LSF系数）和par->lsf（当前帧的LSF系数），
		// 结合插值因子intfact和LPC的阶数LPC_ORD，来生成新的LSF系数数组lsf。
		// 插值操作在语音编码，在处理帧间过渡时，平滑地过渡不同帧之间的LSF系数，从而避免音质上的突变。
		lpc_lsp2pred(lsf, &(lpc[1]), LPC_ORD);
		// 这行代码调用lpc_lsp2pred函数，将插值后的LSF系数lsf转换为对应的LPC滤波器系数，
		// 并将这些系数存储在lpc数组的第二个元素开始的位置（通常lpc[0]是不使用的，或者用于存储其他信息）。
		/* Check signal probability for adaptive spectral enhancement filter */
		temp1 = melpe_add(noise_gain, X12_Q8);
		temp2 = melpe_add(noise_gain, X30_Q8);
		/* sig_prob in Q15 */
		sig_prob = lin_int_bnd(gain, temp1, temp2, 0, ONE_Q15);

		/* Calculate adaptive spectral enhancement filter coefficients */
		ase_num[0] = ONE_Q12;	/* ase_num and ase_den in Q12 */
		temp1 = melpe_mult(sig_prob, ASE_NUM_BW_Q15);
		lpc_bw_expand(&(lpc[1]), &(ase_num[1]), temp1, LPC_ORD);
		temp1 = melpe_mult(sig_prob, ASE_DEN_BW_Q15);
		lpc_bw_expand(&(lpc[1]), ase_den, temp1, LPC_ORD);

		/*      tilt_cof[1] = sig_prob * (intfact * curr_tilt +
		   (1.0 - intfact) * prev_tilt);               */
		temp1 = melpe_mult(curr_tilt, intfact);
		intfact1 = melpe_sub(ONE_Q15, intfact);
		temp2 = melpe_mult(prev_tilt, intfact1);
		temp1 = melpe_add(temp1, temp2);
		tilt_cof[1] = melpe_mult(sig_prob, temp1);	/* tilt_cof in Q15 */

		/* interpolate pitch and pulse gain */
		/*      syn_gain = SYN_GAIN * (intfact * lpc_gain +
		   (1.0 - intfact) * prev_lpc_gain); */
		temp1 = melpe_mult(lpc_gain, intfact);	/* lpc_gain in Q15 */
		temp2 = melpe_mult(prev_lpc_gain, intfact1);
		temp1 = melpe_add(temp1, temp2);	/* temp1 in Q15 */
		syn_gain = melpe_mult(SYN_GAIN_Q4, temp1);
		/* syn_gain in Q4 */
		/*      pitch = intfact * par->pitch + (1.0 - intfact) * prev_par.pitch;  */
		temp1 = melpe_mult(par->pitch, intfact);
		temp2 = melpe_mult(prev_par.pitch, intfact1);
		pitch = melpe_add(temp1, temp2);	/* pitch in Q7 */

		/*      pulse_gain = syn_gain * sqrt(pitch); */
		temp1 = sqrt_fxp(pitch, 7);
		L_temp1 = melpe_L_mult(syn_gain, temp1);
		L_temp1 = melpe_L_shl(L_temp1, 4);
		pulse_gain = melpe_extract_h(L_temp1);	/* pulse_gain in Q0 */

		/* interpolate pulse and noise coefficients */
		temp1 = sqrt_fxp(ifact, 15);
		interp_array(prev_pcof, curr_pcof, pulse_cof, temp1,
			     MIX_ORD + 1);
		interp_array(prev_ncof, curr_ncof, noise_cof, temp1,
			     MIX_ORD + 1);

		set_fc(prev_par.bpvc, &fc_prev);
		set_fc(par->bpvc, &fc_curr);
		temp2 = melpe_sub(ONE_Q15, temp1);
		temp1 = melpe_mult(temp1, fc_curr);	/* temp1 is now Q3 */
		temp2 = melpe_mult(temp2, fc_prev);	/* Q3 */
		fc = melpe_add(temp1, temp2);	/* Q3 */

		/* interpolate jitter */
		/*      jitter = ifact * par->jitter + (1.0 - ifact) * prev_par.jitter;   */
		temp1 = melpe_mult(par->jitter, ifact);
		temp2 = melpe_sub(ONE_Q15, ifact);	/* temp2 is Q15 */
		temp2 = melpe_mult(prev_par.jitter, temp2);
		jitter = melpe_add(temp1, temp2);	/* jitter is Q15 */

		/* scale gain by 0.05 but keep gain in log. */
		/*      gain = pow(10.0, 0.05 * gain); */
		gain = melpe_mult(X005_Q19, gain);	/* gain in Q12 */

		/* Set period length based on pitch and jitter */
		rand_num(&temp1, ONE_Q15, 1);
		/*      length = pitch * (1.0 - jitter * temp) + 0.5; */
		temp1 = melpe_mult(jitter, temp1);
		temp1 = melpe_shr(temp1, 1);	/* temp1 in Q14 */
		temp1 = melpe_sub(ONE_Q14, temp1);
		temp1 = melpe_mult(pitch, temp1);	/* length in Q6 */
		length = melpe_shift_r(temp1, -6);	/* length in Q0 with r_ounding */
		if (length < PITCHMIN)
			length = PITCHMIN;
		if (length > PITCHMAX)
			length = PITCHMAX;

		fill(fs_real, ONE_Q13, length);
		fs_real[0] = 0;
		interp_array(prev_par.fs_mag, par->fs_mag, &fs_real[1], intfact,
			     NUM_HARM);

		harm_syn_pitch(fs_real, &sigbuf[BEGIN], fc, length);
		v_scale(&sigbuf[BEGIN], pulse_gain, length);	/* sigbuf[] is Q0 */

		/* Adaptive spectral enhancement */
		v_equ(&sigbuf[BEGIN - LPC_ORD], ase_del, LPC_ORD);
		lpc_synthesis(&sigbuf[BEGIN], &sigbuf[BEGIN], ase_den, LPC_ORD,
			      length);
		v_equ(ase_del, &sigbuf[BEGIN + length - LPC_ORD], LPC_ORD);

		zerflt(&sigbuf[BEGIN], ase_num, &sigbuf[BEGIN], LPC_ORD,
		       length);
		v_equ(&sigbuf[BEGIN - TILT_ORD], tilt_del, TILT_ORD);
		v_equ(tilt_del, &sigbuf[length + BEGIN - TILT_ORD], TILT_ORD);
		zerflt_Q(&sigbuf[BEGIN], tilt_cof, &sigbuf[BEGIN], TILT_ORD,
			 length, 15);

		/* Possible Signal overflow at this point! */
		/* Perform LPC synthesis filtering */
		v_equ(&sigbuf[BEGIN - LPC_ORD], lpc_del, LPC_ORD);
		lpc_synthesis(&sigbuf[BEGIN], &sigbuf[BEGIN], &(lpc[1]),
			      LPC_ORD, length);
		v_equ(lpc_del, &sigbuf[length + BEGIN - LPC_ORD], LPC_ORD);
		/* Adjust scaling of synthetic speech, sigbuf in Q0 */
		scale_adj(&sigbuf[BEGIN], gain, length, SCALEOVER,
			  INV_SCALEOVER_Q18);

		/* Implement pulse dispersion filter on output speech */
		v_equ(&sigbuf[BEGIN - DISP_ORD], disp_del, DISP_ORD);
		v_equ(disp_del, &sigbuf[length + BEGIN - DISP_ORD], DISP_ORD);
		zerflt_Q(&sigbuf[BEGIN], disp_cof, &sigbuf[BEGIN], DISP_ORD,
			 length, 15);

		/* Copy processed speech to output array (not past frame end) */
		if (melpe_add(syn_begin, length) >= FRAME) {
			v_equ(&sp_out[syn_begin], &sigbuf[BEGIN],
			      (int16_t) (FRAME - syn_begin));

#if POSTFILTER	// sc1200定义为 TRUE , 执行 postfilt
			postfilt(sp_out, prev_par.lsf, par->lsf);
#endif

			/* past end: save remainder in sigsave[0] */
			v_equ(sigsave, &sigbuf[BEGIN + FRAME - syn_begin],
			      (int16_t) (length - (FRAME - syn_begin)));
		} else
			v_equ(&sp_out[syn_begin], &sigbuf[BEGIN], length);

		/* Update syn_begin for next period */
		syn_begin = melpe_add(syn_begin, length);	// 定义 syn_begin
	}

	/* Save previous pulse and noise filters for next frame */
	/* 将当前帧的参数保存为到前一帧的参数里面，以供下一帧语音合成时使用*/
	v_equ(prev_pcof, curr_pcof, MIX_ORD + 1);	// 
	v_equ(prev_ncof, curr_ncof, MIX_ORD + 1);

	/* Copy current parameters to previous parameters for next time */
	prev_par = *par;
	prev_tilt = curr_tilt;
	prev_lpc_gain = lpc_gain;

	/* Update syn_begin for next frame */
	/*syn_begin 更新为当前合成周期的起始位置，以准备开始合成下一帧的语音*/
	syn_begin = melpe_sub(syn_begin, FRAME);
}

/* =========================================================== */
/* melp_syn_init() performs initialization for melp synthesis. */
/* =========================================================== */
// 之前在想第一帧是没有值的，能用到的值都是从这里来的
void melp_syn_init()	// melp_syn_init	
{
	register int16_t i;
	int16_t temp;

	v_zap(prev_par.gain, NUM_GAINFR);
	prev_par.pitch = UV_PITCH_Q7;
	temp = 0;
	for (i = 0; i < LPC_ORD; i++) {
		temp = melpe_add(temp, INV_LPC_ORD);	// INV_LPC_ORD 为 2979
		prev_par.lsf[i] = temp;		// 累加常量
	}
	prev_par.jitter = 0;
	v_zap(&prev_par.bpvc[0], NUM_BANDS);
	syn_begin = 0;	// 在初始化里被定义
	v_zap(sigsave, PITCHMAX);

	fill(prev_par.fs_mag, ONE_Q13, NUM_HARM);

	/* Initialize fixed MSE weighting and inverse of weighting */

	if (!w_fs_init) {
		vq_fsw(w_fs, NUM_HARM, X60_Q9);
		for (i = 0; i < NUM_HARM; i++)
			w_fs_inv[i] = melpe_divide_s(ONE_Q13, w_fs[i]);	/* w_fs_inv in Q14 */	// 被初始化了的，只定义一次
		w_fs_init = TRUE;
	}
}
