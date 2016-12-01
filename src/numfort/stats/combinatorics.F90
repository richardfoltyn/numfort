module numfort_stats_combinatorics

    use, intrinsic :: iso_fortran_env
    use random, only: random_order

    implicit none
    private

contains

subroutine random_order_int32 (x)
    integer, parameter :: INTSIZE = int32
    include "include/random_order_impl.f90"
end subroutine

subroutine random_order_int64 (x)
    integer, parameter :: INTSIZE = int64
    include "include/random_order_impl.f90"
end subroutine

end module
