
#include "numfort.h"

module numfort_optimize_newton

    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_common_input_checks
    use numfort_core, only: signum
    use numfort_optimize_result
    implicit none
    private

    public :: root_newton, root_halley

    integer, parameter :: MSG_LEN = 100

    abstract interface
        subroutine fcn_real32 (x, fx, fpx, args)
            import real32
            integer, parameter :: PREC = real32
            real (PREC), intent(in) :: x
            real (PREC), intent(out) :: fx, fpx
            real (PREC), intent(in out), dimension(:), optional :: args
        end subroutine

        subroutine fcn_real64 (x, fx, fpx, args)
            import real64
            integer, parameter :: PREC = real64
            real (PREC), intent(in) :: x
            real (PREC), intent(out) :: fx, fpx
            real (PREC), intent(in out), dimension(:), optional :: args
        end subroutine

        subroutine fcn_der2_real32 (x, fx, fpx, fppx, args)
            import real32
            integer, parameter :: PREC = real32
            real (PREC), intent(in) :: x
            real (PREC), intent(out) :: fx, fpx, fppx
            real (PREC), intent(in out), dimension(:), optional :: args
        end subroutine

        subroutine fcn_der2_real64 (x, fx, fpx, fppx, args)
            import real64
            integer, parameter :: PREC = real64
            real (PREC), intent(in) :: x
            real (PREC), intent(out) :: fx, fpx, fppx
            real (PREC), intent(in out), dimension(:), optional :: args
        end subroutine

    end interface

    interface root_newton
        module procedure root_newton_real32, root_newton_real64
    end interface

    interface root_halley
        module procedure root_halley_real32, root_halley_real64
    end interface

    interface root_newton_impl
        module procedure root_newton_impl_real32, root_newton_impl_real64
    end interface

    interface check_inputs
        module procedure check_inputs_real32, check_inputs_real64
    end interface

    interface set_defaults
        module procedure set_defaults_real32, set_defaults_real64
    end interface

contains

#define __PREC real32
#include "root_newton_impl.F90"
#undef __PREC

#define __PREC real64
#include "root_newton_impl.F90"
#undef __PREC


end module
