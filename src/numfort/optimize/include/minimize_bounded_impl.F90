

subroutine __APPEND(minimize_bounded,__PREC) (fcn, a, b, x, maxfun, xtol, res)

    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss,__PREC)) :: fcn
        !*  Objective function
    real (PREC), intent(in) :: a
        !*  Lower bound of minimizer domain.
    real (PREC), intent(in) :: b
        !*  Upper bound of minimizer domain.
    real (PREC), intent(out) :: x
        !*  Minimizer on interval [a,b], if one is found. Otherwise, contains
        !   the last guess before terminating.
    integer, intent(in), optional :: maxfun
        !*  Max. number of function evaluations.
    real (PREC), intent(in), optional :: xtol
        !*  Absolute tolerance for termination.
    type (__APPEND(optim_result,__PREC)), intent(in out), optional :: res
        !*  Optimization result.

    type (__APPEND(fwrapper_ss,__PREC)) :: fcn_wrapper
    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res
    integer :: lmaxfun
    real (PREC) :: lxtol

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call check_input (a, b, maxfun, xtol, ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    lxtol = 1e-5_PREC
    lmaxfun = 500
    if (present(xtol)) lxtol = xtol
    if (present(maxfun)) lmaxfun = maxfun

    call wrap_procedure (fcn_wrapper, fcn=fcn)

    call minimize_bounded_impl (fcn_wrapper, a, b, lmaxfun, lxtol, res=ptr_res)
    x = res%x(1)

100 continue
    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

end subroutine


subroutine __APPEND(minimize_bounded_args,__PREC) (fcn, a, b, x, args, maxfun, xtol, res)

    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_args,__PREC)) :: fcn
        !*  Objective function
    real (PREC), intent(in) :: a
        !*  Lower bound of minimizer domain.
    real (PREC), intent(in) :: b
        !*  Upper bound of minimizer domain.
    real (PREC), intent(out) :: x
        !*  Minimizer on interval [a,b], if one is found. Otherwise, contains
        !   the last guess before terminating.
        !*  Array of additional arguments passed to objective function.
    class (args_data), intent(inout) :: args
    integer, intent(in), optional :: maxfun
        !*  Max. number of function evaluations.
    real (PREC), intent(in), optional :: xtol
        !*  Absolute tolerance for termination.
    type (__APPEND(optim_result,__PREC)), intent(in out), optional :: res
        !*  Optimization result.

    type (__APPEND(fwrapper_ss,__PREC)) :: fcn_wrapper
    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res
    integer :: lmaxfun
    real (PREC) :: lxtol

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call check_input (a, b, maxfun, xtol, ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    lxtol = 1e-5_PREC
    lmaxfun = 500
    if (present(xtol)) lxtol = xtol
    if (present(maxfun)) lmaxfun = maxfun

    call wrap_procedure (fcn_wrapper, fcn_args=fcn, args=args)

    call minimize_bounded_impl (fcn_wrapper, a, b, lmaxfun, lxtol, ptr_res)
    x = res%x(1)

100 continue
    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

end subroutine


subroutine __APPEND(check_input,__PREC) (a, b, maxfun, xtol, res)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: a
    real (PREC), intent(in) :: b
    integer, intent(in), optional :: maxfun
    real (PREC), intent(in), optional :: xtol
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res

    res%status = NF_STATUS_OK
    res%msg = ""

    if (b < a) then
        res%status = NF_STATUS_INVALID_ARG
        res%msg = "Invalid bracket, b > a"
        goto 100
    end if

    call check_positive (1, maxfun, "maxfun", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (1.0_PREC, xtol, "xtol", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

100 continue
end subroutine



subroutine __APPEND(minimize_bounded_impl,__PREC) (fcn, xlb, xub, maxfun, &
        xtol, res)
    !*  MINIMIZE_BOUNDED_IMPL minimizes a scalar function on a bounded interval.
    !   Fortran port of Scipy's fminbound/_minimize_scalar_bounded.

    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_ss,__PREC)), intent(in out) :: fcn
    real (PREC), intent(in) :: xlb, xub
    integer, intent(in) :: maxfun
    real (PREC), intent(in) :: xtol
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res

    real (PREC) :: a, b, e, golden_mean, nfc, xf, x, xm, sqrt_eps, rat
    real (PREC) :: fulc, ffulc, fnfc, fx, fu
    real (PREC) :: p, q, r
    real (PREC) :: tol1, tol2
    logical :: golden

    res%status = NF_STATUS_OK
    res%msg = "Solution found"

    sqrt_eps = sqrt(epsilon(0.0_PREC))
    golden_mean = 0.5_PREC * (3.0_PREC - sqrt(5.0_PREC))

    a = xlb
    b = xub

    fulc = a + golden_mean * (b - a)
    nfc = fulc
    xf = fulc
    rat = 0.0_PREC
    e = 0.0_PREC
    x = xf
    call dispatch (fcn, x, fx)

    ffulc = fx
    fnfc = fx
    xm = 0.5_PREC * (a + b)
    tol1 = sqrt_eps * abs(xf) + xtol / 3.0_PREC
    tol2 = 2.0_PREC * tol1

    do while (abs(xf - xm) > (tol2 - 0.5_PREC * (b - a)))

        if (fcn%nfev >= maxfun) then
            res%status = NF_STATUS_MAX_EVAL
            res%status = res%status + NF_STATUS_NOT_CONVERGED
            res%msg = "Max. number of function call reached"
            goto 100
        end if

        golden = .true.
        ! Check for parabolic fit
        if (abs(e) > tol1) then
            golden = .false.
            r = (xf - nfc) * (fx - ffulc)
            q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r
            q = 2.0_PREC * (q - r)
            if (q > 0.0_PREC) p = -p
            q = abs(q)
            r = e
            e = rat

            ! Check for acceptability of parabola
            if ((abs(p) < abs(0.5_PREC*q*r)) .and. (p > q*(a - xf)) .and. &
                    (p < q * (b - xf))) then
                rat = p / q
                x = xf + rat

                if (((x - a) < tol2) .or. ((b - x) < tol2)) then
                    rat = sign(tol1, xm-xf)
                end if
            else
                ! do a golden section step
                golden = .true.
            end if
        end if

        if (golden) then
            ! Do a golden-section step
            if (xf >= xm) then
                e = a - xf
            else
                e = b - xf
            end if
            rat = golden_mean * e
        end if

        x = xf + max(abs(rat), tol1) * sign(1.0_PREC, rat)
        call dispatch (fcn, x, fu)

        if (fu <= fx) then
            if (x >= xf) then
                a = xf
            else
                b = xf
            end if
            fulc = nfc
            ffulc = fnfc
            nfc = xf
            fnfc = fx
            xf = x
            fx = fu
        else
            if (x < xf) then
                a = x
            else
                b = x
            end if

            if ((fu <= fnfc) .or. (nfc == xf)) then
                fulc = nfc
                ffulc = fnfc
                nfc = x
                fnfc = fu
            else if ((fu <= ffulc) .or. (fulc == xf) .or. (fulc == nfc)) then
                fulc = x
                ffulc = fu
            end if
        end if

        xm = 0.5_PREC * (a + b)
        tol1 = sqrt_eps * abs(xf) + xtol / 3.0_PREC
        tol2 = 2.0_PREC * tol1
    end do

100 continue

    call result_update (res, xf, fx, nit=fcn%nfev, nfev=fcn%nfev)

end subroutine
