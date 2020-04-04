      module minpack_real64
        use iso_fortran_env
        implicit none
        private

        integer, parameter :: PREC = real64

        abstract interface
          subroutine lmder_fcn_real64 (m,n,x,fvec,fjac,ldfjac,iflag)
              integer :: m, n, ldfjac, iflag
              double precision :: x(n), fvec(m), fjac(ldfjac, n)
          end subroutine

          subroutine lmdif_fcn_real64 (m, n, x, fvec, iflag)
              integer :: m, n, iflag
              double precision :: x(n), fvec(m)
          end subroutine

          subroutine hybrd_fcn_real64 (n, x, fvec, iflag)
              integer :: n, iflag
              double precision :: x(n), fvec(n)
          end subroutine

          subroutine hybrj_fcn_real64 (n,x,fvec,fjac,ldfjac,iflag)
              integer :: n, iflag, ldfjac
              double precision  :: x(n), fvec(n), fjac(ldfjac, n)
          end subroutine
        end interface

        public :: lmder
        public :: lmdif
        public :: hybrd
        public :: hybrj
        public :: lmder_fcn_real64
        public :: lmdif_fcn_real64
        public :: hybrd_fcn_real64
        public :: hybrj_fcn_real64
        public :: chkder

        real (PREC), parameter :: one   = 1.0_PREC
        real (PREC), parameter :: zero  = 0.0_PREC
        real (PREC), parameter :: p1    = 0.1_PREC
        real (PREC), parameter :: p5    = 0.5_PREC
        real (PREC), parameter :: p25   = 0.25_PREC
        real (PREC), parameter :: p75   = 0.75_PREC
        real (PREC), parameter :: p05   = 0.05_PREC
        real (PREC), parameter :: p001  = 1.0e-3_PREC
        real (PREC), parameter :: p0001 = 1.0e-4_PREC

      contains

      ! driver routines
#include "hybrd.f"
#include "hybrj.f"
#include "lmder.f"
#include "lmdif.f"

      ! helper routines
#include "chkder.f"
#include "dogleg.f"
#include "enorm.f"
#include "fdjac1.f"
#include "fdjac2.f"
#include "lmpar.f"
#include "qform.f"
#include "qrfac.f"
#include "qrsolv.f"
#include "r1mpyq.f"
#include "r1updt.f"



      end module
