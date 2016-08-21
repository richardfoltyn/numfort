module numfort_interpolate_common

    use iso_fortran_env
    implicit none
    private

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
