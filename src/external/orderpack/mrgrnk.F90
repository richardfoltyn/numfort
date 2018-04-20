

module orderpack_mrgrnk
    !*  Module contains (modified) implementation of ORDERPACK's MRGRNK.
    !
    !   Changes applied to the original code:
    !       1.  the implementation was moved to an include file to eliminate
    !           the need to replicate the code one-for-one for each
    !           argument/kind combination.
    !       2.  Additional specific routines for int64 arguments were added.

    use, intrinsic :: iso_fortran_env

    implicit none
    private

    public :: mrgrnk

    interface mrgrnk
        procedure mrgrnk_real32_int32, mrgrnk_real32_int64, &
            mrgrnk_real64_int32, mrgrnk_real64_int64, &
            mrgrnk_int32, mrgrnk_int64
    end interface

    contains

pure subroutine mrgrnk_real32_int32 (XDONT, IRNGT)
    !*  MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.

    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    real (PREC), dimension (:), intent(in) :: XDONT
    integer (INTSIZE), dimension (:), intent(out) :: IRNGT

    real (PREC) :: XVALA, XVALB

#include "mrgrnk_impl.F90"
end subroutine


pure subroutine mrgrnk_real64_int32 (XDONT, IRNGT)
    !*  MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.

    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    real (PREC), dimension (:), intent(in) :: XDONT
    integer (INTSIZE), dimension (:), intent(out) :: IRNGT

    real (PREC) :: XVALA, XVALB

#include "mrgrnk_impl.F90"
end subroutine

pure subroutine mrgrnk_real32_int64 (XDONT, IRNGT)
    !*  MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.

    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int64

    real (PREC), dimension (:), intent(in) :: XDONT
    integer (INTSIZE), dimension (:), intent(out) :: IRNGT

    real (PREC) :: XVALA, XVALB

#include "mrgrnk_impl.F90"
end subroutine


pure subroutine mrgrnk_real64_int64 (XDONT, IRNGT)
    !*  MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.

    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int64

    real (PREC), dimension (:), intent(in) :: XDONT
    integer (INTSIZE), dimension (:), intent(out) :: IRNGT

    real (PREC) :: XVALA, XVALB

#include "mrgrnk_impl.F90"
end subroutine


pure subroutine mrgrnk_int32 (XDONT, IRNGT)
    !*  MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.

    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), dimension (:), intent(in) :: XDONT
    integer (INTSIZE), dimension (:), intent(out) :: IRNGT

    integer (INTSIZE) :: XVALA, XVALB

#include "mrgrnk_impl.F90"
end subroutine


pure subroutine mrgrnk_int64 (XDONT, IRNGT)
    !*  MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.

    integer, parameter :: INTSIZE = int64

    integer (INTSIZE), dimension (:), intent(in) :: XDONT
    integer (INTSIZE), dimension (:), intent(out) :: IRNGT

    integer (INTSIZE) :: XVALA, XVALB

#include "mrgrnk_impl.F90"
end subroutine

end module
