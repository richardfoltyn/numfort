      module lbfgsb_bmnz_real64

        use, intrinsic :: iso_fortran_env, only: real64

        implicit none
        private

        public :: setulb

      contains

#include "lbfgsb.f"
#include "linpack.f"

      end module
