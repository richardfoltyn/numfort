
#include "numfort.h"

module numfort_optimize_newton

    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_optim_result
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

contains

pure subroutine newton_check_inputs (xtol, tol, maxiter, status, msg)
    !*  NEWTON_CHECK_INPUTS performs input validation for both the
    !   plain Newton and the Halley root-finding algorithms.
    !   Used for both real32 and real64 precision.

    integer, parameter :: PREC = real64
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (status_t), intent(out) :: status
    character (*), intent(out) :: msg

    status = NF_STATUS_OK

    if (present(xtol)) then
        if (xtol <= 0.0_PREC) then
            msg = "Tolerance too small"
            goto 100
        end if
    end if

    if (present(tol)) then
        if (tol <= 0.0_PREC) then
            msg = "Tolerance too small"
            goto 100
        end if
    end if

    if (present(maxiter)) then
        if (maxiter <= 0) then
            msg = "maxiter too small"
            goto 100
        end if
    end if

    return

100 continue
    status = NF_STATUS_INVALID_ARG
end subroutine

#define __PREC real32
#include "root_newton_impl.F90"
#undef __PREC

#define __PREC real64
#include "root_newton_impl.F90"
#undef __PREC


end module
