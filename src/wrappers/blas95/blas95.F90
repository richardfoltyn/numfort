

#include <numfort.h>


module blas95
    !*  Module contains Fortran 95 wrappers for BLAS routines that have an
    !   interface which is compatible with the Intel MKL blas95 module.
    !
    !   Note: This module should be used only when Intel MKL is not available,
    !   as otherwise the MKL wrapper can be used directly.

    use, intrinsic :: iso_fortran_env
    use blas_interfaces

    implicit none

    private

    public :: gemv

    public :: gemm


    interface gemv
        procedure gemv_real32, gemv_real64
    end interface

    interface gemm
        procedure gemm_real32, gemm_real64
    end interface

    contains


#include <numfort_real32.h>
#include "gemv.F90"
#include "gemm.F90"


#include <numfort_real64.h>
#include "gemv.F90"
#include "gemm.F90"


end module
