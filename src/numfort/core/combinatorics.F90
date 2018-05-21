
module numfort_core_combinatorics
    !*  Module implements various combinatorial routines such as
    !   the factorial function and n-choose-k.

    use iso_fortran_env, only : real32, real64, int32, int64


    implicit none
    private

    public :: comb
    public :: factorial
    public :: poch

    interface factorial
        module procedure factorial_int32, factorial_int64
    end interface

    interface comb
        module procedure comb_int32, comb_int64
    end interface

    interface poch
        procedure poch_int32, poch_int64
    end interface

contains

!-------------------------------------------------------------------------------
! Factorial function

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



!-------------------------------------------------------------------------------
! COMB routine

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
! Pochhammer symbol (rising factorial)

elemental function poch_int32 (n, k) result(res)
    !*  POCH returns the Pochhammer symbol (rising factorial), defined as
    !   poch(n, k) = n * (n+1) * ... * (n+k-1)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: n
    integer (INTSIZE), intent(in) :: k
    integer (INTSIZE) :: res

    integer (INTSIZE) :: i

    ! For k <= 0, poch(n,k) = 1 (empty product)
    res = 1
    do i = 0, k-1
        res = res * (n + i)
    end do
end function


elemental function poch_int64 (n, k) result(res)
    !*  POCH returns the Pochhammer symbol (rising factorial), defined as
    !   poch(n, k) = n * (n+1) * ... * (n+k-1)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in) :: n
    integer (INTSIZE), intent(in) :: k
    integer (INTSIZE) :: res

    integer (INTSIZE) :: i

    ! For k <= 0, poch(n,k) = 1 (empty product)
    res = 1
    do i = 0, k-1
        res = res * (n + i)
    end do
end function


end module
