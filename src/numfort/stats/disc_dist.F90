#ifdef __INTEL_COMPILER
#define __SUPPORTS_PDT
#endif

#define __DEFAULT_REAL_KIND real64

#ifdef __SUPPORTS_PDT
#define __KIND_PARAM_DECL (PREC)
#define __REAL_KIND_DECL integer, kind :: PREC = __DEFAULT_REAL_KIND
#define __RWRK_KIND PREC
#else
#define __KIND_PARAM_DECL
#define __REAL_KIND_DECL
#define __RWRK_KIND __DEFAULT_REAL_KIND
#endif

module numfort_stats_disc_dist

    use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
    implicit none
    private


    type, abstract :: disc_dist
    contains
        ! pmf method
        procedure, private, pass :: pmf_scalar_i+nt32
        procedure (disc_dist_func_array_int32), deferred, private, pass :: pmf_array_int32
        generic, public :: pmf => pmf_scalar_int32, pmf_array_int32

        ! CDF method
        procedure, private, pass :: cdf_scalar_int32
        procedure (disc_dist_func_array_int32), deferred, private, pass :: cdf_array_int32
        generic, public :: cdf => cdf_scalar_int32, cdf_array_int32

        ! Random sample generation method
        procedure, private, pass :: rvs_scalar_int32
        procedure (disc_dist_array_int32), deferred, private, pass :: rvs_array_int32
        generic, public :: rvs => rvs_scalar_int32, rvs_array_int32

    end type

    interface
        subroutine disc_dist_func_array_int32 (self, x, fx)
            import int32
            import disc_dist
            integer, parameter :: PREC = real64
            integer, parameter :: INTSIZE = int32
            class (disc_dist), intent(in) :: self
            include "include/disc_dist_func_array_args.f90"
        end subroutine

        subroutine disc_dist_array_int32 (self, x)
            import int32
            import disc_dist
            integer, parameter :: INTSIZE = int32
            class (disc_dist), intent(in) :: self
            include "include/disc_dist_array_args.f90"
        end subroutine
    end interface

    public :: disc_dist

contains


! ------------------------------------------------------------------------------
! pmf Method implementation

subroutine pmf_scalar_real32 (self, x, fx)
    integer, parameter :: PREC = real32
    class (disc_dist), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%pmf (x1, fx1)
    fx = fx1(1)
end subroutine

subroutine pmf_scalar_real64 (self, x, fx)
    integer, parameter :: PREC = real64
    class (disc_dist), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%pmf (x1, fx1)
    fx = fx1(1)
end subroutine

! ------------------------------------------------------------------------------
! CDF method implementation

subroutine cdf_scalar_real32 (self, x, fx)
    integer, parameter :: PREC = real32
    class (disc_dist), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%cdf (x1, fx1)
    fx = fx1(1)
end subroutine

subroutine cdf_scalar_real64 (self, x, fx)
    integer, parameter :: PREC = real64
    class (disc_dist), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%cdf (x1, fx1)
    fx = fx1(1)
end subroutine

! ------------------------------------------------------------------------------
! RVS method implementation

subroutine rvs_scalar_real32 (self, x)
    integer, parameter :: PREC = real32
    class (disc_dist), intent(in) :: self
    real (PREC), intent(out) :: x

    real (PREC), dimension(1) :: x1
    x1(1) = x
    call self%rvs (x1)
    x = x1(1)
end subroutine

subroutine rvs_scalar_real64 (self, x)
    integer, parameter :: PREC = real64
    class (disc_dist), intent(in) :: self
    real (PREC), intent(out) :: x

    real (PREC), dimension(1) :: x1
    x1(1) = x
    call self%rvs (x1)
    x = x1(1)
end subroutine

end module
