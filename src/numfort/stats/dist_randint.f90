
#include "dist_pdt.h"

module numfort_stats_dist_randint

    use iso_fortran_env
    use numfort_stats_disc_dist

    implicit none
    private


    type, extends(disc_dist) :: dist_randint
        integer (__PDT_INT_KIND) :: low = 0
        integer (__PDT_INT_KIND) :: high = 1
    contains
        procedure, private, pass :: pmf_array_int32_real64
        procedure, private, pass :: pmf_array_int64_real64

        procedure, private, pass :: rvs_array_int32_real64
        procedure, private, pass :: rvs_array_int64_real64
        procedure, private, pass :: rvs_array_int32_real64_args
        procedure, private, pass :: rvs_array_int64_real64_args
        generic, public :: rvs => rvs_array_int32_real64_args, &
            rvs_array_int64_real64_args
    end type

    type (dist_randint), parameter :: randint = dist_randint(low=0, high=1)

    public :: randint, dist_randint

contains

! ------------------------------------------------------------------------------
! PDF method

subroutine pmf_array_int32_real64 (self, x, fx)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32
    class (dist_randint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: fx

end subroutine

subroutine pmf_array_int64_real64 (self, x, fx)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int64
    class (dist_randint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: fx

end subroutine

! ------------------------------------------------------------------------------
! RVS method
subroutine rvs_array_int32_real64(self, x)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32
    class (dist_randint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(out), dimension(:) :: x

end subroutine

subroutine rvs_array_int64_real64(self, x)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int64
    class (dist_randint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(out), dimension(:) :: x

end subroutine

subroutine rvs_array_int32_real64_args(self, x, low, high)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32
    class (dist_randint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(out), dimension(:) :: x
    integer (INTSIZE), intent(in) :: low
    integer (INTSIZE), intent(in) :: high

end subroutine

subroutine rvs_array_int64_real64_args(self, x, low, high)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int64
    class (dist_randint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(out), dimension(:) :: x
    integer (INTSIZE), intent(in) :: low
    integer (INTSIZE), intent(in) :: high

end subroutine

end module
