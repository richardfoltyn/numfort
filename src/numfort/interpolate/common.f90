module numfort_interpolate_common

    use iso_fortran_env
    use numfort_common, only: ENUM_KIND
    implicit none
    private

    integer (ENUM_KIND), parameter :: INTERP_EVAL_EXTRAPOLATE = 0
    integer (ENUM_KIND), parameter :: INTERP_EVAL_ZERO = 1
    integer (ENUM_KIND), parameter :: INTERP_EVAL_ERROR = 2
    integer (ENUM_KIND), parameter :: INTERP_EVAL_BOUNDARY = 3

    public :: INTERP_EVAL_EXTRAPOLATE, INTERP_EVAL_ZERO, &
        INTERP_EVAL_ERROR, INTERP_EVAL_BOUNDARY

    interface bsearch
        module procedure bsearch_real64, bsearch_real32
    end interface

    interface interp_find
        module procedure interp_find_real32, interp_find_real64
    end interface

    public :: interp_find, bsearch

    contains

pure function bsearch_real64 (needle, haystack) result(lb)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    include "include/bsearch_impl.f90"
end function

pure function bsearch_real32 (needle, haystack) result(lb)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    include "include/bsearch_impl.f90"
end function

pure function interp_find_real64 (needle, haystack) result (res)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    include "include/interp_find_impl.f90"

end function

pure function interp_find_real32 (needle, haystack) result (res)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    include "include/interp_find_impl.f90"

end function



end module
