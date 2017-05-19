
/*
    Include file for real64 precision. 
    THIS FILE IS GENERATED AUTOMATICALLY BY gen_prec_headers.py, DO NOT MODIFY.
*/


#ifdef __PREC
#undef __PREC
#endif
#define __PREC real64

/* BLAS routines */

#ifdef __GEMM
#undef __GEMM
#endif
#define __GEMM DGEMM

/* LAPACK routines */

#ifdef __GESVD
#undef __GESVD
#endif
#define __GESVD DGESVD


#ifdef __GELSD
#undef __GELSD
#endif
#define __GELSD DGELSD


#ifdef __GESDD
#undef __GESDD
#endif
#define __GESDD DGESDD

