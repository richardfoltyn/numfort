module numfort_core

    use iso_fortran_env, only : real32, real64, int32, int64


    implicit none
    private

    public :: comb, signum
    public :: factorial
    public :: PI, PI_real32, PI_real64

    interface factorial
        module procedure factorial_int32, factorial_int64
    end interface

    interface comb
        module procedure comb_int32, comb_int64
    end interface

    interface signum
        module procedure signum_real32, signum_real64, signum_int32, signum_int64
    end interface

    ! Math constants
    real (real64), parameter :: PI_real64 = 3.141592653589793238462643383279502884d0
    real (real32), parameter :: PI_real32 = 3.141592653589793238462643383279502884

    real (real64), parameter :: PI = PI_real64

contains

elemental function factorial_int32(n) result(res)
    integer, parameter :: PREC = int32
    integer (PREC), intent(in) :: n
    integer (PREC) :: res
    integer (PREC) :: i

    res = 1
    do i = 2, n
        res = res * i
    end do
end function

elemental function factorial_int64(n) result(res)
    integer, parameter :: PREC = int64
    integer (PREC), intent(in) :: n
    integer (PREC) :: res
    integer (PREC) :: i

    res = 1
    do i = 2, n
        res = res * i
    end do
end function

elemental function comb_int32 (n, k, repetition) result(res)
    !*  COMB returns the number of combinations (with or without repetition)
    !   of N items taken K at a time.
    integer, parameter :: PREC = int32
#include "comb_impl.F90"
end function

elemental function comb_int64 (n, k, repetition) result(res)
    !*  COMB returns the number of combinations (with or without repetition)
    !   of N items taken K at a time.
    integer, parameter :: PREC = int64
#include "comb_impl.F90"
end function


!-------------------------------------------------------------------------------
! SIGNUM: alternative SIGN function
! Implement behavior that corresponds to what is found in numpy.

elemental function signum_real64 (x) result(res)
    !*  SIGNUM returns -1, 0 or 1 depending on whether x < 0, x == 0 or x > 0.
    !   Note: Differs from intrinsic SIGN, as SIGN return 1 whenever x >= 0
    !   and -1 otherwise.
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: x
    real (PREC) :: res

    if (x > 0.0_PREC) then
        res = 1.0_PREC
    else if (x < 0.0_PREC) then
        res = -1.0_PREC
    else
        res = 0.0_PREC
    end if
end function


elemental function signum_real32 (x) result(res)
    !*  SIGNUM returns -1, 0 or 1 depending on whether x < 0, x == 0 or x > 0.
    !   Note: Differs from intrinsic SIGN, as SIGN return 1 whenever x >= 0
    !   and -1 otherwise.
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: x
    real (PREC) :: res

    if (x > 0.0_PREC) then
        res = 1.0_PREC
    else if (x < 0.0_PREC) then
        res = -1.0_PREC
    else
        res = 0.0_PREC
    end if
end function


elemental function signum_int32 (x) result(res)
    !*  SIGNUM returns -1, 0 or 1 depending on whether x < 0, x == 0 or x > 0.
    !   Note: Differs from intrinsic SIGN, as SIGN return 1 whenever x >= 0
    !   and -1 otherwise.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: x
    integer (INTSIZE)  :: res

    if (x > 0_INTSIZE) then
        res = 1_INTSIZE
    else if (x < 0_INTSIZE) then
        res = -1_INTSIZE
    else
        res = 0_INTSIZE
    end if
end function


elemental function signum_int64 (x) result(res)
    !*  SIGNUM returns -1, 0 or 1 depending on whether x < 0, x == 0 or x > 0.
    !   Note: Differs from intrinsic SIGN, as SIGN return 1 whenever x >= 0
    !   and -1 otherwise.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in) :: x
    integer (INTSIZE)  :: res

    if (x > 0_INTSIZE) then
        res = 1_INTSIZE
    else if (x < 0_INTSIZE) then
        res = -1_INTSIZE
    else
        res = 0_INTSIZE
    end if
end function



end module
