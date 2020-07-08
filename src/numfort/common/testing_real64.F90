

module numfort_common_testing_real64
    !*  Module NUMFORT_COMMON_TESTING contains routines to perform various
    !   tests, eg. if two values are approximately close.

    use, intrinsic :: iso_fortran_env

    use numfort_common_shape

    implicit none

    private

    public :: all_close
    public :: all_close_fast
    public :: is_close

    integer, parameter :: PREC = real64

    interface all_close
        procedure all_close_0d, all_close_1d, all_close_2d, all_close_3d
    end interface

    interface is_close
        procedure is_close
    end interface

    interface all_close_fast
        procedure all_close_fast_1d, all_close_fast_2d, all_close_fast_3d
    end interface

    contains

#include "include/testing_impl.F90"

end module
