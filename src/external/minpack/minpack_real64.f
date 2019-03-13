      module minpack_real64
        use iso_fortran_env
        implicit none
        private

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
