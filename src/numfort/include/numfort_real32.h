/*
    Define macros that transform code to real64 KIND
*/
#ifdef __PREC
#undef __PREC
#undef __APPEND_PREC
#undef __IDENTITY
#endif

#define __PREC real32
#define __IDENTITY(x) x
#define __APPEND_PREC(name) __IDENTITY(name)_real32
