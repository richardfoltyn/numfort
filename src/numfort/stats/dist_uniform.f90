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

module numfort_stats_uniform

    use, intrinsic :: iso_fortran_env
    use numfort_stats_cont_dist, only: cont_dist

    implicit none
    private

    !>  Implementation of the continuous univariate uniform distribution
    !   on an interval [a,b]
    type, extends(cont_dist) :: duniform
        real (real64) :: loc = 0.0d0
            !!  Location parameter such that a = loc
        real (real64) :: scale = 1.0d0
            !!  Scale parameter such that b = loc + scale
    contains
        procedure, private, pass :: pdf_array_real32
        procedure, private, pass :: pdf_array_real64

        procedure, private, pass :: cdf_array_real32
        procedure, private, pass :: cdf_array_real64

        procedure, private, pass :: rvs_array_real32
        procedure, private, pass :: rvs_array_real64

        procedure, private, pass :: get_params_real32
        procedure, private, pass :: get_params_real64
        generic :: get_params => get_params_real32, get_params_real64
    end type

    public :: uniform

    !>  Instance of the standard uniform distribution on [0,1].
    type (duniform), protected :: uniform = duniform(loc=0.0d0, scale=1.0d0)

contains

! ------------------------------------------------------------------------------
! GET_PARAMS method

subroutine get_params_real64 (self, loc, scale, loc_out, scale_out)
    integer, parameter :: PREC = real64
    class (duniform), intent(in) :: self
    real (PREC), intent(in), optional :: loc, scale
    real (PREC), intent(out) :: loc_out, scale_out

    if (present(loc)) then
        loc_out = loc
    else
        loc_out = real(self%loc, PREC)
    end if

    if (present(scale)) then
        scale_out = scale
    else
        scale_out = real(self%scale, PREC)
    end if
end subroutine

subroutine get_params_real32 (self, loc, scale, loc_out, scale_out)
    integer, parameter :: PREC = real32
    class (duniform), intent(in) :: self
    real (PREC), intent(in), optional :: loc, scale
    real (PREC), intent(out) :: loc_out, scale_out

    if (present(loc)) then
        loc_out = loc
    else
        loc_out = real(self%loc, PREC)
    end if

    if (present(scale)) then
        scale_out = scale
    else
        scale_out = real(self%scale, PREC)
    end if
end subroutine

! ------------------------------------------------------------------------------
! PDF method
subroutine pdf_array_real32 (self, x, fx, loc, scale, shape)
    integer, parameter :: PREC = real32
    class (duniform), intent(in) :: self
    include "include/cont_dist_func_array_args.f90"

    real (PREC) :: lloc, lscale

    call self%get_params (loc, scale, lloc, lscale)
    fx = 1 / lscale
    ! clip PDF value outside of valid range
    where (x < lloc)
        fx = 0.0_PREC
    else where (x > (lloc + lscale))
        fx = 0.0_PREC
    end where
end subroutine


subroutine pdf_array_real64 (self, x, fx, loc, scale, shape)
    integer, parameter :: PREC = real64
    class (duniform), intent(in) :: self
    include "include/cont_dist_func_array_args.f90"

    real (PREC) :: lloc, lscale

    call self%get_params (loc, scale, lloc, lscale)
    fx = 1 / lscale
    ! clip PDF value outside of valid range
    where (x < lloc)
        fx = 0.0_PREC
    else where (x > (lloc + lscale))
        fx = 0.0_PREC
    end where
end subroutine

! ------------------------------------------------------------------------------
! CDF method
subroutine cdf_array_real32 (self, x, fx, loc, scale, shape)
    integer, parameter :: PREC = real32
    class (duniform), intent(in) :: self
    include "include/cont_dist_func_array_args.f90"

    real (PREC) :: lloc, lscale

    call self%get_params (loc, scale, lloc, lscale)

    fx = (x - lloc) / lscale
    ! adjust CDF values outside of valid range
    where (x < lloc)
        fx = 0.0_PREC
    else where (x > (lloc + lscale))
        fx = 1.0_PREC
    end where
end subroutine

subroutine cdf_array_real64 (self, x, fx, loc, scale, shape)
    integer, parameter :: PREC = real64
    class (duniform), intent(in) :: self
    include "include/cont_dist_func_array_args.f90"

    real (PREC) :: lloc, lscale

    call self%get_params (loc, scale, lloc, lscale)

    fx = (x - lloc) / lscale
    ! adjust CDF values outside of valid range
    where (x < lloc)
        fx = 0.0_PREC
    else where (x > (lloc + lscale))
        fx = 1.0_PREC
    end where

end subroutine

! ------------------------------------------------------------------------------
! RVS method

subroutine rvs_array_real32 (self, x, loc, scale, shape)
    integer, parameter :: PREC = real32
    class (duniform), intent(in) :: self
    include "include/cont_dist_array_args.f90"

    real (PREC) :: lloc, lscale
    integer :: i

    call self%get_params (loc, scale, lloc, lscale)

    call random_number (x)
    forall (i=1:size(x)) x(i) = x(i) * lscale + lloc
end subroutine

subroutine rvs_array_real64 (self, x, loc, scale, shape)
    integer, parameter :: PREC = real64
    class (duniform), intent(in) :: self
    include "include/cont_dist_array_args.f90"

    real (PREC) :: lloc, lscale
    integer :: i

    call self%get_params (loc, scale, lloc, lscale)

    call random_number (x)
    forall (i=1:size(x)) x(i) = x(i) * lscale + lloc
end subroutine
end module
