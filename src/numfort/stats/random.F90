
module numfort_stats_random
    !*  Module containing procedures for generating random numbers that do not
    !   are not implement in any distribution-specific module.

    use, intrinsic :: iso_fortran_env
    use random, only: random_order_orig => random_order

    implicit none
    private

    public :: random_order
    public :: shuffle
    public :: random_integer

    interface random_order
        procedure random_order_int32, random_order_int64
    end interface

    interface shuffle
        procedure shuffle_int32, shuffle_int64, shuffle_real32, shuffle_real64
    end interface

    interface random_integer
        procedure random_integer_int32, random_integer_int64
    end interface
contains

subroutine random_order_int32 (x, low)
    integer, parameter :: INTSIZE = int32
#include "random_order_impl.f90"
end subroutine

subroutine random_order_int64 (x, low)
    integer, parameter :: INTSIZE = int64
#include "random_order_impl.f90"
end subroutine



subroutine shuffle_int32 (x)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), dimension(:), intent(inout) :: x
    integer (INTSIZE) :: xi
#include "shuffle_impl.F90"
end subroutine

subroutine shuffle_int64 (x)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), dimension(:), intent(inout) :: x
    integer (INTSIZE) :: xi
#include "shuffle_impl.F90"
end subroutine

subroutine shuffle_real32 (x)
    integer, parameter :: PREC = real32
    real (PREC), dimension(:), intent(inout) :: x
    real (PREC) :: xi
#include "shuffle_impl.F90"
end subroutine

subroutine shuffle_real64 (x)
    integer, parameter :: PREC = real64
    real (PREC), dimension(:), intent(inout) :: x
    real (PREC) :: xi
#include "shuffle_impl.F90"
end subroutine



subroutine random_integer_int32 (x, low, high)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:) :: x
    integer (INTSIZE), intent(in) :: low
    integer (INTSIZE), intent(in) :: high
#include "random_integer.F90"
end subroutine

subroutine random_integer_int64 (x, low, high)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:) :: x
    integer (INTSIZE), intent(in) :: low
    integer (INTSIZE), intent(in) :: high
#include "random_integer.F90"
end subroutine

end module
