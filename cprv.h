/* vim: set tabstop=4:softtabstop=4:shiftwidth=4:noexpandtab */

/*------------------------------------------------------------------*/
/*																	*/
/* File:		cprv.h				即 classify.c 头文件！！		*/
/*			但也有别的，比如当前脚本声明结构体在别的函数里用到			*/
/*																	*/
/* Description: Head file for classification, pitch and voicing		*/
/*		仅有一个头文件，用于声明结构体、函数、宏定义！！！				*/
/*-----------------------具体实现由其他源文件完成---------------------*/
/*		比如没有	classify.h	 头文件，想使用就只能include "cprv.h"	*/

#ifndef _CPRV_H_
#define _CPRV_H_

/* ================================================================ *
 *					Definitions										*
 * ================================================================ */
#define EN_UP_RATE_Q15		29491	/* 0.9 * (1 << 15) */
#define	TRACK_NUM			9
#define	CUR_TRACK			2

/* ================================================================ *
 *					Structures										*
 * ================================================================ */

/* ======== Pitch estimation structures ======== */
typedef struct {
	int16_t pit[NODE];	/* integer pitch for each node, Q7 */
	int16_t weight[NODE];	/* time domain correlation, Q15 */
	int16_t cost[NODE];	/* cost function, Q0 */
} pitTrackParam;

typedef struct {
	int16_t classy;	/* the class */
	int16_t subEnergy;	/* full band energy, Q11 */
	int16_t zeroCrosRate;	/* zero crossing rate, Q15 */
	int16_t peakiness;	/* peakiness measure, Q11 */
	int16_t corx;		/* autocorrelation, Q15 */
	int16_t pitch;	/* pitch period, Q0 */
} classParam;

/* ============================ */
/* Prototypes from "classify.c"	*/
/* ============================ */
void classify(int16_t inbuf[], classParam * classStat, int16_t autocorr[]);

#endif
