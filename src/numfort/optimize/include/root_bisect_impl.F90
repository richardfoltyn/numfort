

pure subroutine __APPEND(bisect_check_inputs,__PREC) (a, b, xtol, tol, maxiter, res)
    !*  BISECT_CHECK_INPUTS performs input validation for the bisection
    !   root-finding algorithm.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: a
    real (PREC), intent(in) :: b
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(inout) :: res

    res%status = NF_STATUS_OK

    if (a > b) then
        res%msg = "Invalid bracket, routine requires a < b"
        goto 100
    end if

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




subroutine __APPEND(root_bisect,__PREC) (fcn, a, b, x0, xtol, tol, maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss,__PREC)) :: fcn
    real (PREC), intent(in) :: a
    real (PREC), intent(in) :: b
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(out), optional :: x0
    type (__APPEND(optim_result,__PREC)), intent(out), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn)

    call root_bisect_impl (fwrapper, a, b, x0, xtol, tol, maxiter, res)

end subroutine



subroutine __APPEND(root_bisect_args,__PREC) (fcn, a, b, args, x0, xtol, tol, &
        maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_args,__PREC)) :: fcn
    real (PREC), intent(in) :: a
    real (PREC), intent(in) :: b
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(out), optional :: x0
    type (__APPEND(optim_result,__PREC)), intent(out), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_args=fcn, args=args)

    call root_bisect_impl (fwrapper, a, b, x0, xtol, tol, maxiter, res)

end subroutine



subroutine __APPEND(root_bisect_impl,__PREC) (fcn, a, b, x0, xtol, tol, &
        maxiter, res)
    !*  ROOT_BISECT_IMPL implements a bisection root-finding algorithm.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_ss,__PREC)), intent(inout) :: fcn
    real (PREC), intent(in) :: a
        !*  Lower bound of bracketing interval
    real (PREC), intent(in) :: b
        !*  Upper bound of bracketing interval
    real (PREC), intent(in), optional :: xtol
        !*  Termination criterion in terms of changes in X values in function domain.
        !   Routine will exit when abs(X_n - X_{n-1}) < XTOL.
    real (PREC), intent(in), optional :: tol
        !*  Tolerance level. Routine will exit when abs(f(X)) < TOL
    integer, intent(in), optional :: maxiter
        !*  Maximum number of iterations
    real (PREC), intent(out), optional :: x0
        !*  If present, contains the root on successful exit.
    type (__APPEND(optim_result,__PREC)), intent(out), optional :: res
        !*  Result object

    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res
    real (PREC) :: lxtol, ltol
    integer :: lmaxiter

    real (PREC) :: sgn, x, fx, slb, sub, xlb, xub, flb, fub, xlast
    integer :: iter

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    ! Don't leave uninitialized as these will be passed to RESULT_UPDATE
    ! even if nothing is computed.
    ptr_res%msg = ''
    iter = 0
    x = 0.0
    fx = huge(0.0_PREC)

    call bisect_check_inputs (a, b, xtol, tol, maxiter, res=ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    lmaxiter = 50
    ltol = 1.0e-6_PREC
    lxtol = 1.0e-8_PREC
    if (present(tol)) ltol = tol
    if (present(xtol)) lxtol = xtol
    if (present(maxiter)) lmaxiter = maxiter

    xlb = a
    xub = b

    call dispatch (fcn, xlb, flb)
    call dispatch (fcn, xub, fub)

    slb = signum (flb)
    sub = signum (fub)

    if (slb == sub) then
        ptr_res%msg = 'Initial bounds do not bracket root'
        ptr_res%status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    xlast = xlb

    do iter = 0, lmaxiter - 1
        x = (xub + xlb) / 2.0_PREC
        call dispatch (fcn, x, fx)

        if (abs(fx) < ltol) then
            ptr_res%msg = 'Convergence achieved: abs(f(x)) < TOL'
            ptr_res%status = NF_STATUS_OK
            goto 100
        end if

        if (abs(x-xlast) < lxtol) then
            ptr_res%msg = 'Convergence achieved: abs. change in x small than XTOL'
            ptr_res%status = NF_STATUS_OK
            goto 100
        end if

        ! Update bounds
        sgn = signum (fx)
        if (sgn == slb) then
            xlb = x
            flb = fx
        else
            xub = x
            fub = fx
        end if

        xlast = x
    end do

    ptr_res%msg = 'Max. number of iterations exceeded'
    ptr_res%status = NF_STATUS_MAX_ITER

100 continue

    call result_update (ptr_res, x, fx, nit=iter, nfev=fcn%nfev)

    if (present(x0)) x0 = x

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

end subroutine
