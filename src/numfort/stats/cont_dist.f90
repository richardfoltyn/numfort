module numfort_stats_cont_dist

    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    private


    type, abstract :: cont_dist
    contains
        ! PDF method
        procedure, private, pass :: pdf_scalar_real32
        procedure, private, pass :: pdf_scalar_real64
        procedure (cont_dist_func_array_real32), deferred, private, pass :: pdf_array_real32
        procedure (cont_dist_func_array_real64), deferred, private, pass :: pdf_array_real64
        generic, public :: pdf => pdf_scalar_real32, pdf_scalar_real64, &
            pdf_array_real32, pdf_array_real64

        ! CDF method
        procedure, private, pass :: cdf_scalar_real32
        procedure, private, pass :: cdf_scalar_real64
        procedure (cont_dist_func_array_real32), deferred, private, pass :: cdf_array_real32
        procedure (cont_dist_func_array_real64), deferred, private, pass :: cdf_array_real64
        generic, public :: cdf => cdf_scalar_real32, cdf_scalar_real64, &
            cdf_array_real32, cdf_array_real64

        ! Random sample generation method
        procedure, private, pass :: rvs_scalar_real32
        procedure, private, pass :: rvs_scalar_real64
        procedure (cont_dist_array_real32), deferred, private, pass :: rvs_array_real32
        procedure (cont_dist_array_real64), deferred, private, pass :: rvs_array_real64
        generic, public :: rvs => rvs_scalar_real32, rvs_scalar_real64, &
            rvs_array_real32, rvs_array_real64

    end type

    interface
        subroutine cont_dist_func_array_real32 (self, x, fx, loc, scale, shape)
            import real32
            import cont_dist
            integer, parameter :: PREC = real32
            class (cont_dist), intent(in) :: self
            include "include/cont_dist_func_array_args.f90"
        end subroutine

        subroutine cont_dist_func_array_real64 (self, x, fx, loc, scale, shape)
            import real64
            import cont_dist
            integer, parameter :: PREC = real64
            class (cont_dist), intent(in) :: self
            include "include/cont_dist_func_array_args.f90"
        end subroutine

        subroutine cont_dist_array_real32 (self, x, loc, scale, shape)
            import real32
            import cont_dist
            integer, parameter :: PREC = real32
            class (cont_dist), intent(in) :: self
            include "include/cont_dist_array_args.f90"
        end subroutine

        subroutine cont_dist_array_real64 (self, x, loc, scale, shape)
            import real64
            import cont_dist
            integer, parameter :: PREC = real64
            class (cont_dist), intent(in) :: self
            include "include/cont_dist_array_args.f90"
        end subroutine
    end interface

    public :: cont_dist

contains


! ------------------------------------------------------------------------------
! PDF Method implementation

subroutine pdf_scalar_real32 (self, x, fx, loc, scale, shape)
    integer, parameter :: PREC = real32
    class (cont_dist), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(in), optional :: loc, scale, shape

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%pdf (x1, fx1, loc, scale, shape)
    fx = fx1(1)
end subroutine

subroutine pdf_scalar_real64 (self, x, fx, loc, scale, shape)
    integer, parameter :: PREC = real64
    class (cont_dist), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(in), optional :: loc, scale, shape

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%pdf (x1, fx1, loc, scale, shape)
    fx = fx1(1)
end subroutine

! ------------------------------------------------------------------------------
! CDF method implementation

subroutine cdf_scalar_real32 (self, x, fx, loc, scale, shape)
    integer, parameter :: PREC = real32
    class (cont_dist), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(in), optional :: loc, scale, shape

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%cdf (x1, fx1, loc, scale, shape)
    fx = fx1(1)
end subroutine

subroutine cdf_scalar_real64 (self, x, fx, loc, scale, shape)
    integer, parameter :: PREC = real64
    class (cont_dist), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(in), optional :: loc, scale, shape

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%cdf (x1, fx1, loc, scale, shape)
    fx = fx1(1)
end subroutine

! ------------------------------------------------------------------------------
! RVS method implementation

subroutine rvs_scalar_real32 (self, x, loc, scale, shape)
    integer, parameter :: PREC = real32
    class (cont_dist), intent(in) :: self
    real (PREC), intent(out) :: x
    real (PREC), intent(in), optional :: loc, scale, shape

    real (PREC), dimension(1) :: x1
    x1(1) = x
    call self%rvs (x1, loc, scale, shape)
    x = x1(1)
end subroutine

subroutine rvs_scalar_real64 (self, x, loc, scale, shape)
    integer, parameter :: PREC = real64
    class (cont_dist), intent(in) :: self
    real (PREC), intent(out) :: x
    real (PREC), intent(in), optional :: loc, scale, shape

    real (PREC), dimension(1) :: x1
    x1(1) = x
    call self%rvs (x1, loc, scale, shape)
    x = x1(1)
end subroutine

end module
