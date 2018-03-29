

#include <numfort.h>


module numfort_common_workspace

    use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64

    implicit none
    private

    public :: assert_alloc
    public :: assert_alloc_ptr
    public :: assert_dealloc_ptr
    public :: workspace_get_ptr
    public :: workspace_reset
    public :: workspace_finalize

#include <numfort_real32.h>
#include "workspace_spec.F90"

#include <numfort_real64.h>
#include "workspace_spec.F90"

    contains

#include <numfort_real32.h>
#include "workspace_impl.F90"

#include <numfort_real64.h>
#include "workspace_impl.F90"


end module
