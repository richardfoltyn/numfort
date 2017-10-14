
#include <numfort.h>

module numfort_common_testing
    !*  Module NUMFORT_COMMON_TESTING contains routines to perform various
    !   tests, eg. if two values are approximately close.

    use, intrinsic :: iso_fortran_env
    
    implicit none
    
    private
    
    public :: all_close
    
    interface all_close
        module procedure all_close_real32, all_close_real64, &
            all_close_1d_real32, all_close_1d_real64
    end interface
    
    real (real64), parameter :: DEFAULT_RTOL_real64 = 10.0 ** (log10(epsilon(0.0_real64))/2)
    real (real32), parameter :: DEFAULT_RTOL_real32 = 10.0 ** (log10(epsilon(0.0_real32))/2)
        !*  Default value for the 'rtol' argument

    contains

#include <numfort_real32.h>
#include "include/testing_impl.F90"

#include <numfort_real64.h>
#include "include/testing_impl.F90"


end module
