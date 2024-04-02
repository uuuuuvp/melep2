/* vim: set tabstop=4:softtabstop=4:shiftwidth=4:noexpandtab */

/*
                  
2.4 kbps MELP Proposed Federal Standard speech coder

Fixed-point C code, version 1.0

Copyright (c) 1998, Texas Instruments, Inc.

*/

/*       spbstd.h   SPB standard header file.      */

#ifndef _MACRO_H_
#define _MACRO_H_

#include <math.h>          // 引用形式区别，math.h 为标准库    所以使用<>
#include "ophmconsts.h"    // 引用用户自定义库的时候，使用   ""

/* OSTYPE-dependent definitions/macros. */

#ifdef SunOS4

/* some standard C function definitions missing from SunOS4 */
extern int fclose(FILE * stream);
extern int fprintf(FILE * stream, const char *format, ...);
extern size_t fread(void *ptr, size_t size, size_t nobj, FILE * stream);
extern int fseek(FILE * stream, long offset, int origin);
extern size_t fwrite(const void *ptr, size_t size, size_t nobj, FILE * stream);
extern int printf(const char *format, ...);
extern long random(void);
extern int sscanf(char *s, const char *format, ...);
extern void rewind(FILE * stream);

#else

#endif

/* ====== */
/* Macros */
/* ====== */

#define Max(a, b)		(((a) > (b)) ? (a) : (b))
#define Min(a, b)		(((a) > (b)) ? (b) : (a))

/* =============== */
/* Debugging Modes */
/* =============== */

void inc_saturation();
#if OVERFLOW_CHECK
#define save_saturation()		temp_saturation = saturation
#define restore_saturation()	saturation = temp_saturation
#else
#define save_saturation()
#define restore_saturation()
#endif

#endif
