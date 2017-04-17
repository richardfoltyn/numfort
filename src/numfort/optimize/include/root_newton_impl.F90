

subroutine __APPEND(root_newton,__PREC) (fcn, x, args, xtol, tol, maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fcn,__PREC)) :: fcn
    real (PREC), intent(in out) :: x
    real (PREC), intent(in out), dimension(:), optional :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(in out), optional :: res

    call root_newton_impl (x, args, xtol, tol, maxiter, res, fcn=fcn)
end subroutine

subroutine __APPEND(root_halley,__PREC) (fcn, x, args, xtol, tol, maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fcn_der2,__PREC)) :: fcn
    real (PREC), intent(in out) :: x
    real (PREC), intent(in out), dimension(:), optional :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(in out), optional :: res

    call root_newton_impl (x, args, xtol, tol, maxiter, res, fcn2=fcn)
end subroutine

subroutine __APPEND(root_newton_impl,__PREC) (x, args, xtol, tol, maxiter, res, fcn, fcn2)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in out) :: x
    real (PREC), intent(in out), dimension(:), optional :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(in out), optional :: res
    procedure (__APPEND(fcn,__PREC)), optional :: fcn
        !!  Objective function that computes f(x) and f'(x) in a single call.
    procedure (__APPEND(fcn_der2,__PREC)), optional :: fcn2
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
    
    if (present(tol)) ltol = tol
    if (present(maxiter)) lmaxiter = maxiter
    if (present(xtol)) lxtol = xtol

    ! Note: there is only one common input validation routine which accepts
    ! double precision arguments
    call newton_check_inputs (real(lxtol, real64), real(ltol, real64), lmaxiter, lstatus, msg)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    do_halley = present(fcn2)

    x0 = x
    fppx = 0.0_PREC
    do iter = 1, lmaxiter
        if (do_halley) then
            call fcn2 (x0, fx, fpx, fppx, args)
        else
            call fcn (x0, fx, fpx, args)
        end if

        if (abs(fx) < ltol) then
            msg = "Convergence achieved; abs(f(x)) < tol"
            lstatus = NF_STATUS_OK
            goto 100
        end if

        if (fpx == 0.0_PREC) then
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
            if (discr < 0.0_PREC) then
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
        call result_update (res, x, fx, lstatus, nit=iter, nfev=iter, msg=msg)
    end if
end subroutine
