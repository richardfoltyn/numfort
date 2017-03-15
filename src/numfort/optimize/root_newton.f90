module numfort_optimize_newton

    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_optim_result_mod
    implicit none
    private

    public :: root_newton, root_halley

    integer, parameter :: PREC = real64
    integer, parameter :: MSG_LEN = 100

    abstract interface
        subroutine fobj_real64 (x, fx, fpx, args)
            import real64
            integer, parameter :: PREC = real64
            real (PREC), intent(in) :: x
            real (PREC), intent(out) :: fx, fpx
            real (PREC), intent(in out), dimension(:), optional :: args
        end subroutine

        subroutine fobj_der2_real64 (x, fx, fpx, fppx, args)
            import real64
            integer, parameter :: PREC = real64
            real (PREC), intent(in) :: x
            real (PREC), intent(out) :: fx, fpx, fppx
            real (PREC), intent(in out), dimension(:), optional :: args
        end subroutine
    end interface

    interface root_newton
        module procedure newton_real64
    end interface

    interface root_halley
        module procedure halley_real64
    end interface

contains

subroutine newton_real64 (fcn, x, args, xtol, tol, maxiter, res)
    procedure (fobj_real64) :: fcn
    real (PREC), intent(in out) :: x
    real (PREC), intent(in out), dimension(:), optional :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (optim_result), intent(in out), optional :: res

    call newton_impl_real64 (x, args, xtol, tol, maxiter, res, fcn)
end subroutine

subroutine halley_real64 (fcn, x, args, xtol, tol, maxiter, res)
    procedure (fobj_der2_real64) :: fcn
    real (PREC), intent(in out) :: x
    real (PREC), intent(in out), dimension(:), optional :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (optim_result), intent(in out), optional :: res

    call newton_impl_real64 (x, args, xtol, tol, maxiter, res, fcn2=fcn)
end subroutine

subroutine newton_impl_real64 (x, args, xtol, tol, maxiter, res, fcn, fcn2)
    real (PREC), intent(in out) :: x
    real (PREC), intent(in out), dimension(:), optional :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (optim_result), intent(in out), optional :: res
    procedure (fobj_real64), optional :: fcn
        !!  Objective function that computes f(x) and f'(x) in a single call.
    procedure (fobj_der2_real64), optional :: fcn2
        !!  Objective function that computes f(x), f'(x) and f''(x) in a single call.

    real (PREC) :: ltol, lxtol, fx, fpx, fppx, x0, discr
    integer :: lmaxiter, iter
    type (status_t) :: lstatus
    character (MSG_LEN) :: msg
    logical :: do_halley

    ! Use scipy defaults
    ltol = 1.48d-8
    lxtol = 1.48d-8
    lmaxiter = 50

    msg = ""
    lstatus = NF_STATUS_OK

    call newton_check_inputs (xtol, tol, maxiter, lstatus, msg)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    if (present(tol)) ltol = tol
    if (present(maxiter)) lmaxiter = maxiter
    if (present(xtol)) lxtol = xtol

    do_halley = present(fcn2)

    x0 = x
    fppx = 0.0_PREC
    do iter = 1, lmaxiter
        if (do_halley) then
            call fcn2 (x0, fx, fpx, fppx)
        else
            call fcn (x0, fx, fpx)
        end if

        if (abs(fx) < ltol) then
            msg = "Convergence achieved; abs(f(x)) < tol"
            lstatus = NF_STATUS_OK
            goto 100
        end if

        if (fpx == 0) then
            msg = "Derivative evaluted to 0"
            lstatus = NF_STATUS_OK
            goto 100
        end if

        if (.not. do_halley) then
            ! Newton step
            x = x0 - fx / fpx
        else
            ! Parabolic Halley's method
            discr = fpx ** 2 - 2 * fx * fppx
            if (discr < 0) then
                x = x0 - fpx / fppx
            else
                x = x0 - 2*fx / (fpx + sign(1.0_PREC, fpx) * sqrt(discr))
            end if
        end if

        ! Exit if tolerance level achieved
        if (abs(x - x0) < lxtol) then
            lstatus = NF_STATUS_OK
            msg = "Convergence achieved: abs(x(n)-x(n-1)) < xtol"
            goto 100
        end if

        ! otherwise update and go to next iteration
        x0 = x
    end do

   msg = "Max. number of iterations exceeded"
   lstatus = NF_STATUS_MAX_ITER
   lstatus = lstatus + NF_STATUS_NOT_CONVERGED

100 continue
    if (present(res)) then
        call res%update (x, fx, lstatus, nit=iter, nfev=iter, msg=msg)
    end if
end subroutine


pure subroutine newton_check_inputs (xtol, tol, maxiter, status, msg)
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

end module
