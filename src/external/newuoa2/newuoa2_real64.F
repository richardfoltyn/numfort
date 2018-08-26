      module newuoa2_real64

        use, intrinsic :: iso_fortran_env

        private
        public :: newuoa2

        integer, parameter :: PREC = real64

        interface newuoa2
            procedure newuoa_h
        end interface

#include "constants.f"
#include "interfaces.f"

        contains

#include "bigden.f"
#include "biglag.f"
#include "newuoa-h.f"
#include "newuob-h.f"
#include "trsapp-h.f"
#include "update.f"

      end module