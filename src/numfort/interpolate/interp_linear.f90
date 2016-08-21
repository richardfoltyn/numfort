module numfort_interpolate_linear

    use iso_fortran_env
    use numfort_interpolate_common
    implicit none
    private

    interface interp
        module procedure interp_real32, interp_real64
    end interface

    public :: interp

contains

pure subroutine interp_real32 (x, xp, fp, fx, ext, left, right)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    include "include/interp_impl.f90"
end subroutine

pure subroutine interp_real64 (x, xp, fp, fx, ext, left, right)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    include "include/interp_impl.f90"
end subroutine

end module
