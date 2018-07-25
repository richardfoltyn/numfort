!>  Collection of routines for linear interpolation.

#include <numfort.h>

module numfort_interpolate_linear

    use, intrinsic :: iso_fortran_env
    use numfort_common_kinds
    use numfort_common_status
    use numfort_interpolate_common
    use numfort_interpolate_search

    implicit none
    private

    public :: interp_linear
    public :: interp_linear_impl
    public :: interp_bilinear

#include <numfort_real32.h>
#include "interp_linear_spec.F90"

#include <numfort_real64.h>
#include "interp_linear_spec.F90"

    contains

#include <numfort_real32.h>
#include "interp_linear_impl.F90"

#include <numfort_real64.h>
#include "interp_linear_impl.F90"


end module
