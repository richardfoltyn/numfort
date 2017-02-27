      module minpack_real64
        use iso_fortran_env
        implicit none
        private

        interface
          subroutine lmder_fcn_real64 (m,n,x,fvec,fjac,ldfjac,iflag)
              integer, intent(in) :: m
              integer, intent(in) :: n
              integer, intent(in) :: ldfjac
              double precision, intent(in) :: x(n)
              double precision, intent(in out) :: fvec(m)
              double precision, intent(in out) :: fjac(ldfjac, n)
              integer, intent(in out) :: iflag
          end subroutine

          subroutine lmdif_fcn_real64 (m, n, x, fvec, iflag)
              integer, intent(in) :: m
              integer, intent(in) :: n
              double precision, intent(in) :: x(n)
              double precision, intent(in out) :: fvec(m)
              integer, intent(in out) :: iflag
          end subroutine

          subroutine hybrd_fcn_real64 (n, x, fvec, iflag)
              integer, intent(in) :: n
              double precision, intent(in) :: x(n)
              double precision, intent(in out) :: fvec(n)
              integer, intent(in out) :: iflag
          end subroutine

          subroutine hybrj_fcn_real64 (n,x,fvec,fjac,ldfjac,iflag)
              integer, intent(in) :: n
              integer, intent(in) :: ldfjac
              double precision, intent(in) :: x(n)
              double precision, intent(in out) :: fvec(n)
              double precision, intent(in out) :: fjac(ldfjac, n)
              integer, intent(in out) :: iflag
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


      contains

      ! driver routines
      include "hybrd.f"
      include "hybrj.f"
      include "lmder.f"
      include "lmdif.f"

      ! helper routines
      include "chkder.f"
      include "dogleg.f"
      include "dpmpar.f"
      include "enorm.f"
      include "fdjac1.f"
      include "fdjac2.f"
      include "lmpar.f"
      include "qform.f"
      include "qrfac.f"
      include "qrsolv.f"
      include "r1mpyq.f"
      include "r1updt.f"



      end module
