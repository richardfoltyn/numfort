

module numfort_common_testing_real32
    !*  Module NUMFORT_COMMON_TESTING contains routines to perform various
    !   tests, eg. if two values are approximately close.

    use, intrinsic :: iso_fortran_env

    use numfort_common_shape

    implicit none

    private

    public :: all_close
    public :: is_close

    integer, parameter :: PREC = real32

    interface all_close
        procedure all_close_0d, all_close_1d, all_close_2d, all_close_3d
    end interface

    interface is_close
        procedure is_close
    end interface

    contains

#include "include/testing_impl.F90"

end module
