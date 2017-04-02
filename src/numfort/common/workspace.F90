#include "numfort.h"


module numfort_common_workspace

    use, intrinsic :: iso_fortran_env, only : real32, real64, int32

    implicit none
    private

    public :: workspace_real32, workspace_real64
    public :: assert_alloc, assert_alloc_ptr, assert_dealloc_ptr

    integer, parameter :: SIZE_UNALLOCATED = -1

    type, abstract :: workspace_base
        integer, dimension(:), allocatable :: iwrk
        logical, dimension(:), allocatable :: lwrk
        character (len=:), allocatable :: cwrk
        integer :: nrwrk = SIZE_UNALLOCATED, niwrk = SIZE_UNALLOCATED, &
            ncwrk = SIZE_UNALLOCATED, nlwrk = SIZE_UNALLOCATED
    end type

    type, extends(workspace_base) :: workspace_real64
        real (real64), dimension(:), allocatable :: rwrk
    end type

    type, extends(workspace_base) :: workspace_real32
        real (real32), dimension(:), allocatable :: rwrk
    end type

    interface assert_alloc
        module procedure ws_assert_alloc_real32, ws_assert_alloc_real64
    end interface

    interface assert_alloc_ptr
        module procedure ws_assert_alloc_ptr_real32, ws_assert_alloc_ptr_real64
    end interface

    interface assert_dealloc_ptr
        module procedure ws_assert_dealloc_ptr_real32, ws_assert_dealloc_ptr_real64
    end interface

contains

#define __PREC real32
#include "workspace_impl.F90"
#undef __PREC

#define __PREC real64
#include "workspace_impl.F90"
#undef __PREC


end module
