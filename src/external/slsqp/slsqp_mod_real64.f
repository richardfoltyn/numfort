      module slsqp_mod_real64
        use, intrinsic :: iso_fortran_env

        implicit none
        private

        public :: slsqp
        public :: slsqp_data, linmin_data

        integer, parameter :: PREC = real64

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
