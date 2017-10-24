

pure subroutine __APPEND(check_inputs,__PREC) (xtol, tol, maxiter, res)
    !*  CHECK_INPUTS performs input validation for both the
    !   plain Newton and the Halley root-finding algorithms.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res

    res%status = NF_STATUS_OK

    call check_positive (0.0_PREC, xtol, "xtol", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, tol, "tol", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (0, maxiter, "maxiter", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    return

100 continue
    res%status = NF_STATUS_INVALID_ARG
end subroutine


pure subroutine __APPEND(set_defaults,__PREC) (xtol_in, tol_in, maxiter_in, &
        xtol, tol, maxiter)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), optional :: xtol_in
    real (PREC), intent(in), optional :: tol_in
    integer, intent(in), optional :: maxiter_in
    real (PREC), intent(out) :: xtol, tol
    integer, intent(out) :: maxiter

    ! Use scipy defaults
    tol = sqrt(epsilon(0.0_PREC))
    xtol = sqrt(epsilon(0.0_PREC))
    maxiter = 50

    if (present(tol_in)) tol = tol_in
    if (present(maxiter_in)) maxiter = maxiter_in
    if (present(xtol_in)) xtol = xtol_in

end subroutine


subroutine __APPEND(root_newton,__PREC) (fcn, x, args, xtol, tol, maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fcn,__PREC)) :: fcn
    real (PREC), intent(in out) :: x
    real (PREC), intent(in out), dimension(:), optional :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(in out), optional :: res

    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res

    real (PREC) :: lxtol, ltol
    integer :: lmaxiter

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call check_inputs (xtol, tol, maxiter, ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    call set_defaults (xtol, tol, maxiter, lxtol, ltol, lmaxiter)

    call root_newton_impl (x, args, lxtol, ltol, lmaxiter, ptr_res, fcn=fcn)

100 continue
    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)
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

    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res

    real (PREC) :: lxtol, ltol
    integer :: lmaxiter

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call check_inputs (xtol, tol, maxiter, ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    call set_defaults (xtol, tol, maxiter, lxtol, ltol, lmaxiter)

    call root_newton_impl (x, args, lxtol, ltol, lmaxiter, ptr_res, fcn2=fcn)

100 continue
    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)
end subroutine


subroutine __APPEND(root_newton_impl,__PREC) (x, args, xtol, tol, maxiter, &
        res, fcn, fcn2)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in out) :: x
    real (PREC), intent(in out), dimension(:), optional :: args
    real (PREC), intent(in) :: xtol
    real (PREC), intent(in) :: tol
    integer, intent(in) :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res
    procedure (__APPEND(fcn,__PREC)), optional :: fcn
        !!  Objective function that computes f(x) and f'(x) in a single call.
    procedure (__APPEND(fcn_der2,__PREC)), optional :: fcn2
        !!  Objective function that computes f(x), f'(x) and f''(x) in a single call.

    real (PREC) :: fx, fpx, fppx, x0, discr
    integer :: iter
    logical :: do_halley

    do_halley = present(fcn2)

    x0 = x
    fppx = 0.0_PREC
    do iter = 1, maxiter
        if (do_halley) then
            call fcn2 (x0, fx, fpx, fppx, args)
        else
            call fcn (x0, fx, fpx, args)
        end if

        if (abs(fx) < tol) then
            res%msg = "Convergence achieved; abs(f(x)) < tol"
            res%status = NF_STATUS_OK
            goto 100
        end if

        if (fpx == 0.0_PREC) then
            res%msg = "Derivative evaluted to 0"
            res%status = NF_STATUS_OK
            goto 100
        end if

        if (.not. do_halley) then
            ! Newton step
            x = x0 - fx / fpx
        else
            ! Parabolic Halley's method
            discr = fpx ** 2.0_PREC - 2.0_PREC * fx * fppx
            if (discr < 0.0_PREC) then
                x = x0 - fpx / fppx
            else
                x = x0 - 2.0_PREC*fx / (fpx + signum(fpx) * sqrt(discr))
            end if
        end if

        ! Exit if tolerance level achieved
        if (abs(x - x0) < xtol) then
            res%status = NF_STATUS_OK
            res%msg = "Convergence achieved: abs(x(n)-x(n-1)) < xtol"
            goto 100
        end if

        ! otherwise update and go to next iteration
        x0 = x
    end do

   res%msg = "Max. number of iterations exceeded"
   res%status = NF_STATUS_MAX_ITER
   res%status = res%status + NF_STATUS_NOT_CONVERGED

100 continue
    call result_update (res, x, fx, nit=iter, nfev=iter)
end subroutine
