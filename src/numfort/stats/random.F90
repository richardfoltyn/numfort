
module numfort_stats_random
    !*  Module containing procedures for generating random numbers that do not
    !   are not implement in any distribution-specific module.

    use, intrinsic :: iso_fortran_env
    use random, only: random_order_orig => random_order

    implicit none
    private

    public :: random_order
    public :: random_integer

    interface random_order
        procedure random_order_int32, random_order_int64
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
