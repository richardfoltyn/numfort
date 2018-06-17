

module numfort_interpolate_search

    use, intrinsic :: iso_fortran_env

    implicit none

    private

    public :: interp_find, bsearch


    interface bsearch
        module procedure bsearch_real64, bsearch_real32, bsearch_int32
    end interface

    interface interp_find
        module procedure interp_find_real32, interp_find_real64
    end interface


   contains


pure function bsearch_real64 (needle, haystack) result(lb)
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:) :: haystack

#include "bsearch_impl.F90"
end function


pure function bsearch_real32 (needle, haystack) result(lb)
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:) :: haystack

#include "bsearch_impl.F90"
end function


pure function bsearch_int32 (needle, haystack) result(lb)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: needle
    integer (INTSIZE), intent(in), dimension(:) :: haystack

#include "bsearch_impl.F90"
end function


pure function interp_find_real64 (needle, haystack) result (res)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

#include "interp_find_impl.F90"
end function

pure function interp_find_real32 (needle, haystack) result (res)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

#include "interp_find_impl.F90"
end function



end module