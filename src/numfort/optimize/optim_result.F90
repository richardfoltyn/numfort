! Module implements OPTIM_RESULT, a type used as optimization result
! object for all optimization routines.
! Author: Richard Foltyn

#include "numfort.h"

module numfort_optim_result

    use, intrinsic :: iso_fortran_env, only: real32, real64, int32
    use numfort_common
    use numfort_common_alloc

    implicit none
    private

    public :: optim_result_real32, optim_result_real64
    public :: result_update

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
        module procedure update_scalar_scalar_real32, update_scalar_scalar_real64
    end interface

    interface result_update
        module procedure update_vec_scalar_real32, update_vec_scalar_real64
    end interface

    interface result_update
        module procedure update_real32, update_real64
    end interface

    interface result_reset
        module procedure reset_real32, reset_real64
    end interface

contains

#define __PREC real32
#include "optim_result_impl.F90"
#undef __PREC

#define __PREC real64
#include "optim_result_impl.F90"
#undef __PREC


end module
