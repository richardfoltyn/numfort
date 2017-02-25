      module fitpack_real64

        use iso_fortran_env
        use constraints_tree_mod
        implicit none
        private

        public :: curfit
        public :: concon
        public :: splev
        public :: splder

      contains

      include "concon.f"
      include "curfit.f"
      include "fpback.f"
      include "fpbspl.f"
      include "fpchec.f"
      include "fpcoco.f"
      include "fpcosp.f"
      include "fpcurf.f"
      include "fpdisc.f"
      include "fpgivs.f"
      include "fpknot.f"
      include "fprati.f"
      include "fprota.f"
      include "splev.f"
      include "splder.f"

      end module
