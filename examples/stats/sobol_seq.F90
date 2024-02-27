
#include "numfort.h"

program sobol_demo

    use, intrinsic :: iso_fortran_env
    use numfort_stats, only: sobol_init, sobol_next, sobol_state
    implicit none


    call test_int32_real32 ()
    call test_int32_real64 ()
    call test_int64_real32 ()
    call test_int64_real64 ()

contains

#define __PREC real32

#define __INTSIZE int32
#include "sobol_seq_impl.F90"
#undef __INTSIZE

#define __INTSIZE int64
#include "sobol_seq_impl.F90"
#undef __INTSIZE

#undef __PREC


#define __PREC real64

#define __INTSIZE int32
#include "sobol_seq_impl.F90"
#undef __INTSIZE

#define __INTSIZE int64
#include "sobol_seq_impl.F90"
#undef __INTSIZE

#undef __PREC

end program
