module numfort_interpolate_linear

    use iso_fortran_env
    use numfort_interpolate_common
    implicit none
    private

    interface interp_linear
        module procedure interp_scalar_real32, interp_scalar_real64, &
            interp_array_real32, interp_array_real64
    end interface

    interface interp_linear_impl
        module procedure interp_linear_real32, interp_linear_real64
    end interface

    public :: interp_linear

contains

pure subroutine interp_linear_real32 (x, xp, fp, fx, ext, left, right)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    include "include/interp_linear_impl.f90"
end subroutine

pure subroutine interp_linear_real64 (x, xp, fp, fx, ext, left, right)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    include "include/interp_linear_impl.f90"
end subroutine

pure subroutine interp_scalar_real32 (x, xp, fp, fx, ext, left, right)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    include "include/interp_linear_scalar.f90"
end subroutine

pure subroutine interp_scalar_real64 (x, xp, fp, fx, ext, left, right)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    include "include/interp_linear_scalar.f90"
end subroutine

pure subroutine interp_array_real32 (x, xp, fp, fx, ext, left, right)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    include "include/interp_linear_array.f90"
end subroutine

pure subroutine interp_array_real64 (x, xp, fp, fx, ext, left, right)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    include "include/interp_linear_array.f90"
end subroutine

end module
