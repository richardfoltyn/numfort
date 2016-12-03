/*
    Note: CPP/FPP symbols __XXX__ are definitions added by build system,
    definitions __XXX are defined in the source code.
*/

#ifndef __NUMFORT_DIST_PDT_H__
#define __NUMFORT_DIST_PDT_H__

/* include definitions for building PDTs */
#include <numfort_pdt.h>

/* Real KIND parameter to use for PDTs */
#define __DEFAULT_REAL_KIND real64
/* Default INT paramter to use for PDTs -- use compiler default */
#define __DEFAULT_INT_KIND __NUMFORT_DEFAULT_INT_KIND__

/* To be used when declaring attibutes of PDTs */
#ifdef __NUMFORT_SUPPORTS_PDT_KIND__
#define __PDT_INT_KIND INTSIZE
#define __PDT_REAL_KIND PREC
#else
#define __PDT_INT_KIND __DEFAULT_INT_KIND
#define __PDT_REAL_KIND __DEFAULT_REAL_KIND
#endif

#endif
