#include "numfort.h"


module numfort_common_workspace

    use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64

    implicit none
    private

    public :: workspace_real32, workspace_real64
    public :: assert_alloc, assert_alloc_ptr, assert_dealloc_ptr
    public :: workspace_get_ptr
    public :: workspace_finalize, workspace_reset

    integer, parameter :: SIZE_UNALLOCATED = -1

    type :: workspace_real64
        integer (int64), private :: roffset = 0, ioffset = 0, coffset = 0, loffset = 0
        integer, dimension(:), allocatable :: iwrk
        logical, dimension(:), allocatable :: lwrk
        character (len=:), allocatable :: cwrk
        integer :: nrwrk = SIZE_UNALLOCATED, niwrk = SIZE_UNALLOCATED, &
            ncwrk = SIZE_UNALLOCATED, nlwrk = SIZE_UNALLOCATED
        real (real64), dimension(:), allocatable :: rwrk
    contains
        final :: ws_finalize_real64
    end type

    type :: workspace_real32
        integer (int64), private :: roffset = 0, ioffset = 0, coffset = 0, loffset = 0
        integer, dimension(:), allocatable :: iwrk
        logical, dimension(:), allocatable :: lwrk
        character (len=:), allocatable :: cwrk
        integer :: nrwrk = SIZE_UNALLOCATED, niwrk = SIZE_UNALLOCATED, &
            ncwrk = SIZE_UNALLOCATED, nlwrk = SIZE_UNALLOCATED
        real (real32), dimension(:), allocatable :: rwrk
    contains
        final :: ws_finalize_real32
    end type

    interface assert_alloc
        module procedure ws_assert_alloc_int32_real32, &
            ws_assert_alloc_int32_real64, ws_assert_alloc_int64_real32, &
            ws_assert_alloc_int64_real64
    end interface

    interface assert_alloc_ptr
        module procedure ws_assert_alloc_ptr_real32, ws_assert_alloc_ptr_real64
    end interface

    interface assert_dealloc_ptr
        module procedure ws_assert_dealloc_ptr_real32, ws_assert_dealloc_ptr_real64
    end interface

    interface workspace_finalize
        module procedure ws_finalize_real32, ws_finalize_real64
    end interface

    interface workspace_get_ptr
        module procedure ws_get_ptr_1d_int32_real32, ws_get_ptr_2d_int32_real32, &
            ws_get_ptr_1d_int32_real64, ws_get_ptr_2d_int32_real64
    end interface

    interface workspace_reset
        module procedure ws_reset_real32, ws_reset_real64
    end interface

contains

#include "numfort_real32.h"
#include "workspace_impl.F90"

#include "numfort_real64.h"
#include "workspace_impl.F90"


end module
