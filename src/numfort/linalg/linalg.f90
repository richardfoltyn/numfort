module numfort_linalg

    use iso_fortran_env

    implicit none
    private

    interface inv
        procedure :: inv_real64, inv_real32
    end interface

    public :: inv

contains

subroutine inv_real64(A, Ainv, info)
    integer, parameter :: PREC = real64
    include "include/inv_impl.f90"
end subroutine

subroutine inv_real32(A, Ainv, info)
    integer, parameter :: PREC = real32
    include "include/inv_impl.f90"
end subroutine

end module
