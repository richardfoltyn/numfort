module numfort_interpolate_common

    use iso_fortran_env
    use numfort_common, only: NF_ENUM_KIND
    implicit none
    private

    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_EXTRAPOLATE = 0
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_ZERO = ishft(1, 0)
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_ERROR = ishft(1, 1)
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_BOUNDARY = ishft(1, 2)
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_CONST = ishft(1, 3)

    interface bsearch
        module procedure bsearch_real64, bsearch_real32, bsearch_int32
    end interface

    interface interp_find
        module procedure interp_find_real32, interp_find_real64
    end interface

    public :: interp_find, bsearch

contains

pure function bsearch_real64 (needle, haystack) result(lb)
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:) :: haystack

    include "include/bsearch_impl.f90"
end function


pure function bsearch_real32 (needle, haystack) result(lb)
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:) :: haystack

    include "include/bsearch_impl.f90"
end function


pure function bsearch_int32 (needle, haystack) result(lb)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: needle
    integer (INTSIZE), intent(in), dimension(:) :: haystack

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
