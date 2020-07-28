


module numfort_common_workspace_real32

    use, intrinsic :: iso_fortran_env, only : real32, int32, int64

    use numfort_common_status

    implicit none
    private

    public :: assert_alloc
    public :: assert_alloc_ptr
    public :: assert_dealloc_ptr
    public :: workspace_get_ptr
    public :: workspace_reset

    integer, parameter :: PREC = real32

#include "workspace_spec.F90"

    contains

#include "workspace_impl.F90"


end module
