/* vim: set tabstop=4:softtabstop=4:shiftwidth=4:noexpandtab */

/* ================================================================== */
/*                                                                    */
/*    Microsoft Speech coder     ANSI-C Source Code                   */
/*    SC1200 1200 bps speech coder                                    */
/*    Fixed Point Implementation      Version 7.0                     */
/*    Copyright (C) 2000, Microsoft Corp.                             */
/*    All rights reserved.                                            */
/*                   MELPe 组合前面的脚本写成模块                       */
/* ================================================================== */
/* ========================================= */
/* melp.c: Mixed Excitation LPC speech coder */
/* ========================================= */

//STANAG-4991

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "sc1200.h"
#include "global.h"
#include "macro.h"
#include "constant.h"

#include "mat_lib.h"
#include "mathhalf.h"
#include "dsp_sub.h"		// 语音处理
#include "melp_sub.h"
#include "math_lib.h"
#include "math.h"

#if NPP		// 预处理模式（ 有 # ）根据是否定义了 NPP 的宏，来看是否 包含 npp.h
#include "npp.h"
#endif

#define X05_Q7				64	/* 0.5 * (1 << 7) */
#define THREE_Q7			384	/* 3 * (1 << 7) */

/* ====== External memory ====== */

int16_t mode;
int16_t chwordsize;
int16_t bitBufSize, bitBufSize12, bitBufSize24;
/* ========== Static definations ========== */

#define PROGRAM_NAME			"SC1200 1200 bps speech coder"
#define PROGRAM_VERSION			"Version 7 / 42 Bits"
#define PROGRAM_DATE			"10/25/2000"

/* ========== Public Prototypes ========== */
void melpe_n(short *sp);                       // 降噪
void melpe_i(void);							   // 初始化
void melpe_a(unsigned char *buf, short *sp);   // 压缩
void melpe_s(short *sp, unsigned char *buf);   // 解压缩
void melpe_al(unsigned char *buf, short *sp);		// 2400 ?
void melpe_i2(void); 

//------------------------------NPP------------------------

//denoise 180 samples sp->sp    降噪
void melpe_n(short *sp)
{
	//noise preprocessor for other codecs (frame=180)	2400 
	npp(sp, sp);   // 在melpe_a 里面也有npp！
}

//-------------------------1200----------------------------------

//init melpe codec at 1200 bps    
void melpe_i(void)    
{
	//====== Run MELPE codec ====== 
	mode = ANA_SYN;		// sc1200.h 定义ANA_SYN 为 0 
	rate = RATE2400;
	frameSize = (int16_t) FRAME;	// define BLOCK		(NF*FRAME), 3*180=540

	chwordsize = 8;
	bitNum12 = 81;
	bitNum24 = 54;
	bitBufSize12 = 11;	// 比特数
	bitBufSize24 = 7;
	bitBufSize = bitBufSize24;

	// 前面已说明 MODE 为0
	melp_ana_init();	// 给 analysis 初始化 参数
	melp_syn_init();	// 给 synthesis 初始化 参数
}

//compress 540 samples sp (67.5 mS) -> 81 bits buf (11 bytes)
void melpe_a(unsigned char *buf, short *sp)  // (bitstream,speech)前面的是输出的比特流
{
	npp(sp, sp);
	// npp(&(sp[FRAME]), &(sp[FRAME]));
	// npp(&(sp[2 * FRAME]), &(sp[2 * FRAME]));
	analysis(sp, melp_par);   // wav => speech => melp_par => cnbuf
	memcpy(buf, chbuf, 7);   // cnbuf => buf ( bitstream )
}

//decompress 54 bits buf (7 bytes) -> 180 samples sp (22.5 mS)
void melpe_s(short *sp, unsigned char *buf)  /// 给你 bitstream ，输出speech
{
	//syntesis
	memcpy(chbuf, buf, 7);   //  buf (bitstream) => cnbuf 
	synthesis(melp_par, sp);  //  cnbuf => melp_par => speech		// cnbuf 扮演什么角色？
}