
#include "disc_dist.h"

module numfort_stats_disc_dist

    use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
    implicit none
    private


    type, abstract :: disc_dist __PDT_PARAM_DECL_BOTH (PREC, INTSIZE)
        __PDT_KIND_DECL(PREC, __DEFAULT_REAL_KIND)
        __PDT_KIND_DECL(INTSIZE, __DEFAULT_INT_KIND)
    contains
        ! pmf method
        procedure, private, pass :: pmf_scalar_int32_real64
        procedure, private, pass :: pmf_scalar_int64_real64
        procedure (disc_dist_func_array_int32_real64), deferred, private, pass :: &
            pmf_array_int32_real64
        procedure (disc_dist_func_array_int64_real64), deferred, private, pass :: &
            pmf_array_int64_real64
        generic, public :: pmf => pmf_scalar_int32_real64, pmf_scalar_int64_real64, &
            pmf_array_int32_real64, pmf_array_int64_real64

        ! CDF method
        ! procedure, private, pass :: cdf_scalar
        ! procedure (disc_dist_func_array), deferred, private, pass :: cdf_array
        ! generic, public :: cdf => cdf_scalar, cdf_array

        ! Random sample generation method
        procedure, private, pass :: rvs_scalar_int32_real64
        procedure, private, pass :: rvs_scalar_int64_real64
        procedure (disc_dist_array_int32_real64), deferred, private, pass :: &
            rvs_array_int32_real64
        procedure (disc_dist_array_int64_real64), deferred, private, pass :: &
            rvs_array_int64_real64
        generic, public :: rvs => rvs_scalar_int32_real64, rvs_scalar_int64_real64, &
            rvs_array_int32_real64, rvs_array_int64_real64

    end type

    interface
        subroutine disc_dist_array_int32_real64 (self, x)
            import int32, real64
            import disc_dist
            integer, parameter :: INTSIZE = int32
            integer, parameter :: PREC = real64
            class (disc_dist __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
            integer (INTSIZE), intent(out), dimension(:) :: x
        end subroutine

        subroutine disc_dist_array_int64_real64 (self, x)
            import int64, real64
            import disc_dist
            integer, parameter :: INTSIZE = int64
            integer, parameter :: PREC = real64
            class (disc_dist __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
            integer (INTSIZE), intent(out), dimension(:) :: x
        end subroutine

        subroutine disc_dist_func_array_int32_real64 (self, x, fx)
            import int32, real64
            import disc_dist
            integer, parameter :: PREC = real64
            integer, parameter :: INTSIZE = int32
            class (disc_dist __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
            integer (INTSIZE), intent(in), dimension(:) :: x
            real (PREC), intent(out), dimension(:) :: fx
        end subroutine

        subroutine disc_dist_func_array_int64_real64 (self, x, fx)
            import int64, real64
            import disc_dist
            integer, parameter :: PREC = real64
            integer, parameter :: INTSIZE = int64
            class (disc_dist __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
            integer (INTSIZE), intent(in), dimension(:) :: x
            real (PREC), intent(out), dimension(:) :: fx
        end subroutine
    end interface

    public :: disc_dist

contains


! ------------------------------------------------------------------------------
! pmf Method implementation

subroutine pmf_scalar_int32_real64 (self, x, fx)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32
    class (disc_dist __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(in) :: x
    real (PREC), intent(out) :: fx

    integer (INTSIZE), dimension(1):: x1
    real (PREC), dimension(1) :: fx1
    x1(1) = x
    call self%pmf (x1, fx1)
    fx = fx1(1)
end subroutine

subroutine pmf_scalar_int64_real64 (self, x, fx)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int64
    class (disc_dist __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(in) :: x
    real (PREC), intent(out) :: fx

    integer (INTSIZE), dimension(1) :: x1
    real (PREC), dimension(1) :: fx1
    x1(1) = x
    call self%pmf (x1, fx1)
    fx = fx1(1)
end subroutine

! ------------------------------------------------------------------------------
! CDF method implementation
!
! subroutine cdf_scalar_real32 (self, x, fx)
!     integer, parameter :: PREC = real32
!     class (disc_dist), intent(in) :: self
!     real (PREC), intent(in) :: x
!     real (PREC), intent(out) :: fx
!
!     real (PREC), dimension(1) :: x1, fx1
!     x1(1) = x
!     call self%cdf (x1, fx1)
!     fx = fx1(1)
! end subroutine
!
! subroutine cdf_scalar_real64 (self, x, fx)
!     integer, parameter :: PREC = real64
!     class (disc_dist), intent(in) :: self
!     real (PREC), intent(in) :: x
!     real (PREC), intent(out) :: fx
!
!     real (PREC), dimension(1) :: x1, fx1
!     x1(1) = x
!     call self%cdf (x1, fx1)
!     fx = fx1(1)
! end subroutine
!
! ------------------------------------------------------------------------------
! RVS method implementation

subroutine rvs_scalar_int32_real64 (self, x)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32
    class (disc_dist __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(out) :: x

    integer (INTSIZE), dimension(1) :: x1
    x1(1) = x
    call self%rvs (x1)
    x = x1(1)
end subroutine

subroutine rvs_scalar_int64_real64 (self, x)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: PREC = real64
    class (disc_dist __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
    integer (INTSIZE), intent(out) :: x

    integer (INTSIZE), dimension(1) :: x1
    x1(1) = x
    call self%rvs (x1)
    x = x1(1)
end subroutine

end module
