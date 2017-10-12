      module slsqp_mod_real64
        use, intrinsic :: iso_fortran_env
        
        implicit none
        private
        
        public :: slsqp
        public :: slsqp_data
        
        integer, parameter :: PREC = real64

        ! Constants that are used throughout various subroutines
        real (PREC), parameter :: epmach = epsilon(1.0_PREC)
        real (PREC), parameter :: zero = 0.0_PREC
        real (PREC), parameter :: one = 1.0_PREC
        real (PREC), parameter :: two = 2.0_PREC
        real (PREC), parameter :: four = 4.0_PREC
        real (PREC), parameter :: ten = 10.0_PREC
        real (PREC), parameter :: hun = 100.0_PREC

        type :: slsqp_data
            real (PREC) :: alpha, f0, gs, h1, h2, h3, h4, t, t0, tol
            integer :: iexact, incons, ireset, itermx, line, n1, n2, n3
        end type
        
        contains
        
        include "slsqp_optmz.f"
      end module
