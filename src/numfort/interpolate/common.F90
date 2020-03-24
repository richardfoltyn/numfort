
#include <numfort.h>

module numfort_interpolate_common

    use iso_fortran_env
    use numfort_common, only: NF_ENUM_KIND
    use numfort_common_status

    implicit none
    private

    public :: check_input_ext
    public :: search_cache

    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_EXTRAPOLATE = 0
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_ZERO = ishft(1, 0)
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_ERROR = ishft(1, 1)
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_BOUNDARY = ishft(1, 2)
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_CONST = ishft(1, 3)

    type :: search_cache
        integer :: i = 1
    end type

#include <numfort_real32.h>
#include "common_spec.F90"

#include <numfort_real64.h>
#include "common_spec.F90"

    contains


#include <numfort_real32.h>
#include "common_impl.F90"

#include <numfort_real64.h>
#include "common_impl.F90"


end module
