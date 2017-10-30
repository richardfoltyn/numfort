! Module implements OPTIM_RESULT, a type used as optimization result
! object for all optimization routines.
! Author: Richard Foltyn

#include "numfort.h"

module numfort_optimize_result

    use, intrinsic :: iso_fortran_env, only: real32, real64, int32
    use numfort_arrays_copy
    use numfort_common
    use numfort_common_alloc

    implicit none
    private

    public :: optim_result_real32, optim_result_real64
    public :: result_update, result_reset
    public :: assert_alloc_ptr, assert_dealloc_ptr

    integer, parameter :: UNINITIALIZED_COUNTER = -1

    type :: optim_result_real32
        integer :: nfev = UNINITIALIZED_COUNTER, nit = UNINITIALIZED_COUNTER
        type (status_t) :: status
        logical :: success = .false.
        character (100) :: msg
        real (real32), dimension(:), allocatable :: x, fx
    end type

    type :: optim_result_real64
        integer :: nfev = UNINITIALIZED_COUNTER, nit = UNINITIALIZED_COUNTER
        type (status_t) :: status
        logical :: success = .false.
        character (100) :: msg
        real (real64), dimension(:), allocatable :: x, fx
    end type

    interface result_update
        module procedure update_ss_real32, update_ss_real64
    end interface

    interface result_update
        module procedure update_vs_real32, update_vs_real64
    end interface

    interface result_update
        module procedure update_real32, update_real64
    end interface

    interface result_update
        module procedure update_int_status_real32, update_int_status_real64
    end interface

    interface result_reset
        module procedure reset_real32, reset_real64
    end interface

    interface assert_alloc_ptr
        module procedure assert_alloc_ptr_real32, assert_alloc_ptr_real64
    end interface

    interface assert_dealloc_ptr
        module procedure assert_dealloc_ptr_real32, assert_dealloc_ptr_real64
    end interface

contains

#define __PREC real32
#include "optim_result_impl.F90"
#undef __PREC

#define __PREC real64
#include "optim_result_impl.F90"
#undef __PREC


end module
