      module quadpack_real64
        !*  Wrapper for 64-bit precision version of QUADPACK routines.

        use, intrinsic :: iso_fortran_env

        implicit none

        private

        public :: qags, qagse

        integer, parameter :: PREC = real64

        abstract interface
            function f_integrand (x) result(fx)
                import PREC
                real (PREC), intent(in) :: x
                real (PREC) :: fx
            end function
        end interface

        interface qags
            module procedure dqags
        end interface

        interface qagse
            module procedure dqagse
        end interface

        contains

        include "dqags.f"
        include "dqagse.f"
        include "dqelg.f"
        include "dqk21.f"
        include "dqpsrt.f"

      end module
