/*
    Define macros that transform code to real64 KIND
*/
#ifdef __PREC
#undef __PREC
#undef __APPEND_PREC
#endif

#define __PREC real64
#define __APPEND_PREC(name) name##_real64
