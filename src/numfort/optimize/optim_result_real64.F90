

module numfort_optimize_result_real64
    !*  Module implements OPTIM_RESULT, a type used as optimization result
    !   object for all optimization routines.

    use, intrinsic :: iso_fortran_env, only: real64

    use numfort_common_status
    use numfort_common_copy_alloc, only: copy_alloc

    implicit none

    integer, parameter :: PREC = real64

    private

    public :: result_update
    public :: result_reset
    public :: result_finalize
    public :: assert_alloc_ptr
    public :: assert_dealloc_ptr

    integer, parameter :: UNINITIALIZED_COUNTER = 0

    type, public :: optim_result
        integer :: nfev = UNINITIALIZED_COUNTER
        !*  Number of function evaluations performed
        integer :: nit = UNINITIALIZED_COUNTER
        !*  Number of iterations performed
        type (status_t) :: status
        !*  Detailed exit status
        logical :: success = .false.
        !*  Exit status flag
        character (100) :: msg
        real (PREC), dimension(:), allocatable :: x
        real (PREC), dimension(:), allocatable :: fx
    end type

    interface result_update
        procedure update, update_ss, update_vs
    end interface

    interface result_reset
        procedure reset
    end interface

    interface assert_alloc_ptr
        procedure assert_alloc_ptr
    end interface

    interface assert_dealloc_ptr
        procedure assert_dealloc_ptr
    end interface

    interface result_finalize
        procedure result_finalize
    end interface

    contains

#include "optim_result_impl.F90"

end module
