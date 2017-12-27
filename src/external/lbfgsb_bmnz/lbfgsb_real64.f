      module lbfgsb_bmnz_real64
        use iso_fortran_env
        implicit none
        private

        public :: setulb
      contains
        include "lbfgsb.f"
        include "linpack.f"
      end module
