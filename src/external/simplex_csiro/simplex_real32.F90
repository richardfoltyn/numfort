module simplex_csiro_real32

    use, intrinsic :: iso_fortran_env

    implicit none
    private

    public :: minim
    public :: functn_if

    integer, parameter :: PREC = real32

    abstract interface
        subroutine functn_if (x, fx)
            import PREC
            real (PREC), intent(in), dimension(:), contiguous :: x
            real (PREC), intent(out) :: fx
        end subroutine
    end interface

    contains

#include "minim.f90"

end module
