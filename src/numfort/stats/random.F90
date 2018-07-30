
module numfort_stats_random
    !*  Module containing procedures for generating random numbers that do not
    !   are not implement in any distribution-specific module.

    use, intrinsic :: iso_fortran_env
    use random, only: random_order_orig => random_order

    implicit none
    private

    interface random_order
        module procedure random_order_int32, random_order_int64
    end interface

    public :: random_order

contains

subroutine random_order_int32 (x, low)
    integer, parameter :: INTSIZE = int32
#include "random_order_impl.f90"
end subroutine

subroutine random_order_int64 (x, low)
    integer, parameter :: INTSIZE = int64
#include "random_order_impl.f90"
end subroutine

end module
