
#include <numfort.h>

module numfort_optimize_result
    !*  Module implements OPTIM_RESULT, a type used as optimization result
    !   object for all optimization routines.

    use, intrinsic :: iso_fortran_env, only: real32, real64, int32
    use numfort_common
    use numfort_common_alloc

    implicit none
    private

    public :: result_update
    public :: result_reset
    public :: result_finalize
    public :: assert_alloc_ptr
    public :: assert_dealloc_ptr

    integer, parameter :: UNINITIALIZED_COUNTER = 0

#include <numfort_real32.h>
#include "optim_result_spec.F90"

#include <numfort_real64.h>
#include "optim_result_spec.F90"

    contains

#include <numfort_real32.h>
#include "optim_result_impl.F90"

#include <numfort_real64.h>
#include "optim_result_impl.F90"


end module
