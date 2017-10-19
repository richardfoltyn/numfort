      module slsqp_mod_real64
        use, intrinsic :: iso_fortran_env

        implicit none
        private

        public :: slsqp
        public :: slsqp_data, linmin_data

        integer, parameter :: PREC = real64

C       Interfaces for BLAS 1 routines used in SLSQP implementation
        interface
            function ddot (n, dx, incx, dy, incy)
                import real64
                integer, intent(in) :: n, incx, incy
                real (real64), intent(in) :: dx(*), dy(*)
                real (real64) :: ddot
            end function

            function dnrm2 (n, x, incx)
                import real64
                integer, intent(in) :: incx, n
                real (real64), intent(in) :: x(*)
                real (real64) :: dnrm2
            end function

            subroutine dcopy (n, dx, incx, dy, incy)
                import real64
                integer, intent(in) :: n, incx, incy
                real (real64), intent(in) :: dx(*)
                real (real64), intent(out) :: dy(*)
            end subroutine

            subroutine daxpy (n, da, dx, incx, dy, incy)
                import real64
                integer, intent(in) :: n, incx, incy
                real (real64), intent(in) :: da, dx(*)
                real (real64), intent(in out) :: dy(*)
            end subroutine

            subroutine dscal (n, da, dx, incx)
                import real64
                integer, intent(in) :: n, incx
                real (real64), intent(in) :: da
                real (real64), intent(in out) :: dx(*)
            end subroutine
        end interface


C       Constants that are used throughout various subroutines
        real (PREC), parameter :: epmach = epsilon(1.0_PREC)
        real (PREC), parameter :: zero = 0.0_PREC
        real (PREC), parameter :: one = 1.0_PREC
        real (PREC), parameter :: two = 2.0_PREC
        real (PREC), parameter :: four = 4.0_PREC
        real (PREC), parameter :: ten = 10.0_PREC
        real (PREC), parameter :: hun = 100.0_PREC

        type :: slsqp_data
C           Data container to store variables that were declared using
C           SAVE or initialized in a DATA statement in routine SLSQPB.
            real (PREC) :: alpha    = 0.0_PREC
            real (PREC) :: f0       = 0.0_PREC
            real (PREC) :: gs       = 0.0_PREC
            real (PREC) :: h1       = 0.0_PREC
            real (PREC) :: h2       = 0.0_PREC
            real (PREC) :: h3       = 0.0_PREC
            real (PREC) :: h4       = 0.0_PREC
            real (PREC) :: t        = 0.0_PREC
            real (PREC) :: t0       = 0.0_PREC
            real (PREC) :: tol      = 0.0_PREC
            integer :: iexact       = 0
            integer :: incons       = 0
            integer :: ireset       = 0
            integer :: itermx       = 0
            integer :: line         = 0
            integer :: n1           = 0
            integer :: n2           = 0
            integer :: n3           = 0
        end type

        type :: linmin_data
C           Data container to store variables that were declared using
C           SAVE or initialized in a DATA statement in routine LINMIN.
            real (PREC) :: a  = 0.0_PREC
            real (PREC) :: b  = 0.0_PREC
            real (PREC) :: d  = 0.0_PREC
            real (PREC) :: e  = 0.0_PREC
            real (PREC) :: u  = 0.0_PREC
            real (PREC) :: v  = 0.0_PREC
            real (PREC) :: w  = 0.0_PREC
            real (PREC) :: x  = 0.0_PREC
            real (PREC) :: fu = 0.0_PREC
            real (PREC) :: fv = 0.0_PREC
            real (PREC) :: fw = 0.0_PREC
            real (PREC) :: fx = 0.0_PREC
        end type

        contains

        include "slsqp_optmz.f"
      end module
