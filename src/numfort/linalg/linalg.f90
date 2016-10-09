module numfort_linalg

    use iso_fortran_env
    use lapack_interfaces

    implicit none
    private

    interface inv
        module procedure :: inv_real64, inv_real32
    end interface

    interface det
        module procedure :: det_real32, det_real64
    end interface

    public :: inv, det

contains

subroutine inv_real64(A, Ainv, info)
    integer, parameter :: PREC = real64
    include "include/inv_impl.f90"
end subroutine

subroutine inv_real32(A, Ainv, info)
    integer, parameter :: PREC = real32
    include "include/inv_impl.f90"
end subroutine

function det_real32 (a) result(d)
    integer, parameter :: PREC = real32
    include "include/det_impl.f90"
end function

function det_real64 (a) result(d)
    integer, parameter :: PREC = real64
    include "include/det_impl.f90"
end function

end module
