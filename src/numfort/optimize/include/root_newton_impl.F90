

pure subroutine __APPEND(check_inputs,__PREC) (xtol, tol, maxiter, xtol2, res)
    !*  CHECK_INPUTS performs input validation for both the
    !   plain Newton and the Halley root-finding algorithms.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: xtol2
    type (__APPEND(optim_result,__PREC)), intent(inout) :: res

    res%status = NF_STATUS_OK

    call check_positive (0.0_PREC, xtol, "xtol", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, xtol2, "xtol2", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, tol, "tol", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (0, maxiter, "maxiter", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    return

100 continue
    res%status = NF_STATUS_INVALID_ARG
end subroutine


pure subroutine __APPEND(check_bracket,__PREC) (x, a, b, res)
    !*  CHECK_BRACKET verifies that if A, B are present, A < B and
    !   X >= A, X <= B.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
    real (PREC), intent(in), optional :: a, b
    type (__APPEND(optim_result,__PREC)), intent(inout) :: res

    res%status = NF_STATUS_OK

    if (present(a)) then
        if (a > x) then
            res%msg = 'Invalid bracket lower bound'
            goto 100
        end if
    end if

    if (present(b)) then
        if (b < x) then
            res%msg = 'Invalid bracket upper bound'
            goto 100
        end if
    end if

    if (present(a) .and. present(b)) then
        if (a >= b) then
            res%msg = 'Invalid bracketing interval'
            goto 100
        end if
    end if

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


subroutine __APPEND(root_newton_args,__PREC) (fcn, x, args, ndiff, xtol, tol, &
        maxiter, xtol2, dstep, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_args,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    class (args_data), intent(inout) :: args
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: xtol2
        !*  If present, identify period-2 cycles if |X(n)-X(n-2)| < XTOL2
        !   and take the midpoint of [X(n), X(n-1)] as the next candidate root
        !   whenever such cycles occur.
    real (PREC), intent(in), optional :: dstep
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn_args=fcn, args=args, eps=dstep)

    call root_newton_impl (fwrapper, x, xtol, tol, maxiter, xtol2, res)

end subroutine

subroutine __APPEND(root_newton_jac_args,__PREC) (fcn, jac, x, args, xtol, &
        tol, maxiter, xtol2, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_args,__PREC)) :: fcn
    procedure (__APPEND(fss_args,__PREC)) :: jac
    real (PREC), intent(inout) :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: xtol2
        !*  If present, identify period-2 cycles if |X(n)-X(n-2)| < XTOL2
        !   and take the midpoint of [X(n), X(n-1)] as the next candidate root
        !   whenever such cycles occur.
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_args=fcn, jac_args=jac, args=args)

    call root_newton_impl (fwrapper, x, xtol, tol, maxiter, xtol2, res)

end subroutine

subroutine __APPEND(root_newton_fcn_jac_args,__PREC) (fcn, x, args, xtol, &
        tol, maxiter, xtol2, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_fcn_jac_args,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: xtol2
        !*  If present, identify period-2 cycles if |X(n)-X(n-2)| < XTOL2
        !   and take the midpoint of [X(n), X(n-1)] as the next candidate root
        !   whenever such cycles occur.
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac_args=fcn, args=args)

    call root_newton_impl (fwrapper, x, xtol, tol, maxiter, xtol2, res)

end subroutine


subroutine __APPEND(root_newton,__PREC) (fcn, x, ndiff, xtol, tol, &
        maxiter, xtol2, dstep, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: xtol2
        !*  If present, identify period-2 cycles if |X(n)-X(n-2)| < XTOL2
        !   and take the midpoint of [X(n), X(n-1)] as the next candidate root
        !   whenever such cycles occur.
    real (PREC), intent(in), optional :: dstep
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            call result_reset (res)
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn=fcn, eps=dstep)

    call root_newton_impl (fwrapper, x, xtol, tol, maxiter, xtol2, res)

end subroutine

subroutine __APPEND(root_newton_jac,__PREC) (fcn, jac, x, xtol, &
        tol, maxiter, xtol2, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss,__PREC)) :: fcn
    procedure (__APPEND(fss,__PREC)) :: jac
    real (PREC), intent(inout) :: x
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: xtol2
        !*  If present, identify period-2 cycles if |X(n)-X(n-2)| < XTOL2
        !   and take the midpoint of [X(n), X(n-1)] as the next candidate root
        !   whenever such cycles occur.
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn, jac=jac)

    call root_newton_impl (fwrapper, x, xtol, tol, maxiter, xtol2, res)

end subroutine

subroutine __APPEND(root_newton_fcn_jac,__PREC) (fcn, x, xtol, &
        tol, maxiter, xtol2, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_fcn_jac,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: xtol2
        !*  If present, identify period-2 cycles if |X(n)-X(n-2)| < XTOL2
        !   and take the midpoint of [X(n), X(n-1)] as the next candidate root
        !   whenever such cycles occur.
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac=fcn)

    call root_newton_impl (fwrapper, x, xtol, tol, maxiter, xtol2, res)

end subroutine



subroutine __APPEND(root_newton_impl,__PREC) (fcn, x, xtol, tol, maxiter, &
        xtol2, res)
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_ss,__PREC)), intent(inout) :: fcn
    real (PREC), intent(inout) :: x
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: xtol2
        !*  If present, identify period-2 cycles if |X(n)-X(n-2)| < XTOL2
        !   and take the midpoint of [X(n), X(n-1)] as the next candidate root
        !   whenever such cycles occur.
    type (__APPEND(optim_result,__PREC)), intent(inout), optional, target :: res

    real (PREC) :: fx, fpx, x0
    real (PREC) :: x2
    integer :: iter

    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res
    real (PREC) :: lxtol, ltol
    integer :: lmaxiter
    logical :: do_xtol, check_cycle

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call check_inputs (xtol, tol, maxiter, xtol2, ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    call set_defaults (xtol, tol, maxiter, lxtol, ltol, lmaxiter)

    x0 = x

    x2 = x

    do iter = 1, lmaxiter
        ! Turn ON checking for convergence in terms of XTOL
        do_xtol = .true.

        call dispatch_fcn_jac (fcn, x0, fx, fpx)

        if (abs(fx) < ltol) then
            ptr_res%msg = "Convergence achieved; abs(f(x)) < ltol"
            ptr_res%status = NF_STATUS_OK
            goto 100
        end if

        if (fpx == 0.0_PREC) then
            ptr_res%msg = "Derivative evaluted to 0"
            ptr_res%status = NF_STATUS_OK
            goto 100
        end if

        ! Newton step
        x = x0 - fx / fpx

        ! Check for period-2 cycles
        check_cycle = (modulo(iter-3,2) == 0) .and. iter >= 3 .and. present(xtol2)
        if (check_cycle) then
            if (abs(x-x2) < xtol2) then
                x = (x + x0) / 2.0_PREC
                ! Skip checking convergence in terms of XTOL in this iteration
                do_xtol = .false.
            end if

            ! Store for next cycle check
            x2 = x
        end if

        ! Exit if tolerance level achieved
        if (do_xtol .and. abs(x - x0) < lxtol) then
            ptr_res%status = NF_STATUS_OK
            ptr_res%msg = "Convergence achieved: abs(x(n)-x(n-1)) < lxtol"
            goto 100
        end if

        ! otherwise update and go to next iteration
        x0 = x
    end do

   ptr_res%msg = "Max. number of iterations exceeded"
   ptr_res%status = NF_STATUS_MAX_ITER
   ptr_res%status = ptr_res%status + NF_STATUS_NOT_CONVERGED

100 continue

    call result_update (ptr_res, x, fx, nit=iter, nfev=fcn%nfev)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)
end subroutine


subroutine __APPEND(newton_bisect_args,__PREC) (fcn, x, args, ndiff, a, b, xtol, &
        tol, maxiter, dstep, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_args,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    class (args_data), intent(inout) :: args
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: a, b
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: dstep
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn_args=fcn, args=args, eps=dstep)

    call newton_bisect_impl (fwrapper, x, a, b, xtol, tol, maxiter, res)

end subroutine

subroutine __APPEND(newton_bisect_jac_args,__PREC) (fcn, jac, x, args, a, b, &
        xtol, tol, maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_args,__PREC)) :: fcn
    procedure (__APPEND(fss_args,__PREC)) :: jac
    real (PREC), intent(inout) :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: a, b
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_args=fcn, jac_args=jac, args=args)

    call newton_bisect_impl (fwrapper, x, a, b, xtol, tol, maxiter, res)

end subroutine


subroutine __APPEND(newton_bisect,__PREC) (fcn, x, ndiff, a, b, xtol, tol, &
        maxiter, dstep, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: a, b
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: dstep
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            call result_reset (res)
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn=fcn, eps=dstep)

    call newton_bisect_impl (fwrapper, x, a, b, xtol, tol, maxiter, res)

end subroutine

subroutine __APPEND(newton_bisect_jac,__PREC) (fcn, jac, x, a, b, xtol, &
        tol, maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss,__PREC)) :: fcn
    procedure (__APPEND(fss,__PREC)) :: jac
    real (PREC), intent(inout) :: x
    real (PREC), intent(in), optional :: a, b
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn, jac=jac)

    call newton_bisect_impl (fwrapper, x, a, b, xtol, tol, maxiter, res)

end subroutine


subroutine __APPEND(newton_bisect_fcn_jac_args,__PREC) (fcn, x, a, b, xtol, &
        tol, maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_fcn_jac,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    real (PREC), intent(in), optional :: a, b
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac=fcn)

    call newton_bisect_impl (fwrapper, x, a, b, xtol, tol, maxiter, res)

end subroutine



subroutine __APPEND(newton_bisect_fcn_jac,__PREC) (fcn, x, args, a, b, xtol, &
        tol, maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fss_fcn_jac_args,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    real (PREC), intent(in), optional :: a, b
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_ss,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac_args=fcn, args=args)

    call newton_bisect_impl (fwrapper, x, a, b, xtol, tol, maxiter, res)

end subroutine


subroutine __APPEND(newton_bisect_impl,__PREC) (fcn, x, a, b, xtol, tol, &
        maxiter, res)
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_ss,__PREC)), intent(inout) :: fcn
    real (PREC), intent(inout) :: x
    real (PREC), intent(in), optional :: a, b
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(inout), optional, target :: res

    real (PREC) :: fx, fpx, x0
    real (PREC) :: xlb, xub, s, xstart, slb, sub, dub, dlb
    real (PREC) :: fa, fb
    integer :: iter
    logical :: has_bracket

    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res
    real (PREC) :: lxtol, ltol
    integer :: lmaxiter

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call check_inputs (xtol, tol, maxiter, res=ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    call check_bracket (x, a, b, ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    call set_defaults (xtol, tol, maxiter, lxtol, ltol, lmaxiter)

    xstart = x

    x0 = x
    call dispatch_fcn_jac (fcn, x0, fx, fpx)
    ! Check for immediate convergence before wasting any more computation time
    if (abs(fx) < ltol) then
        ptr_res%msg = "Convergence achieved; abs(f(x)) < ltol"
        ptr_res%status = NF_STATUS_OK
        goto 100
    end if

    ! Create initial bracket. If user provided A, B use those, otherwise
    ! use the starting point X.
    ! If both A and B are missing, then the bracket will start out as a single
    ! point X until we sample more function values. In that case
    ! SLB holds sign of function value until bracket has been identified.
    if (present(a)) then
        call dispatch (fcn, a, fa)
        slb = signum (fa)
        xlb = a
    else
        slb = signum (fx)
        xlb = x
    end if

    if (present(b)) then
        call dispatch (fcn, b, fb)
        sub = signum (fb)
        xub = b
    else
        sub = signum (fx)
        xub = x
    end if

    if (xlb > xub) then
        ! Flip values
        call swap (xlb, xub)
        call swap (slb, sub)
    end if

    ! Bracket is established whenever the signs are different
    has_bracket = slb*sub < 0.0_PREC

    ! If we have no bracket, but A, B have been provided by user, return
    ! error.
    if (.not. has_bracket .and. present(a) .and. present(b)) then
        ptr_res%msg = "Initial interval does not bracket root"
        ptr_res%status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    do iter = 1, lmaxiter

        if (abs(fx) < ltol) then
            ptr_res%msg = "Convergence achieved; abs(f(x)) < ltol"
            ptr_res%status = NF_STATUS_OK
            goto 100
        end if

        if (fpx == 0.0_PREC) then
            ptr_res%msg = "Derivative evaluted to 0"
            ptr_res%status = NF_STATUS_OK
            goto 100
        end if

        ! Newton step
        x = x0 - fx / fpx

        if (has_bracket) then
            if (x < xlb .or. x > xub) then
                x = (xlb + xub) / 2.0_PREC
            end if
        end if

        ! Compute function value and derivative for the NEXT iteration
        call dispatch_fcn_jac (fcn, x, fx, fpx)

        ! Exit if tolerance level achieved.
        ! Note: do this only after updating FX, we need to return correct FX=f(X)
        if (abs(x - x0) < lxtol) then
            ptr_res%status = NF_STATUS_OK
            ptr_res%msg = "Convergence achieved: abs(x(n)-x(n-1)) < lxtol"
            goto 100
        end if

        s = slb * signum (fx)
        if (.not. has_bracket) then
            if (s < 0.0_PREC) then
                ! Create bracket
                if (x > xub) then
                    xlb = xub
                    xub = x
                    ! SLB remains unchanged
                    sub = signum (fx)
                else if (x < xlb) then
                    xub = xlb
                    xlb = x
                    sub = slb
                    slb = signum (fx)
                else
                    dub = abs(xstart - xub)
                    dlb = abs(xstart - xlb)

                    if (dub < dlb) then
                        xlb = x
                        ! xub remains unchanged
                        sub = slb
                        slb = signum (fx)
                    else
                        xub = x
                        sub = signum (fx)
                        ! XLB, SLB remain unchanged
                    end if
                end if

                has_bracket = .true.

            else
                ! Note: if s = 0 then we must have fx = 0 and algorithm will
                ! terminate in the next iteration, so we can ignore that case.

                ! Update boundaries if sign has not changed
                xlb = min(xlb, x)
                xub = max(xub, x)
            end if
        else
            ! Update bracket
            if (s > 0.0_PREC) then
                ! f(x) has same sign as f(xlb)
                xlb = x
            else
                xub = x
            end if
        end if

        x0 = x
    end do

   ptr_res%msg = "Max. number of iterations exceeded"
   ptr_res%status = NF_STATUS_MAX_ITER
   ptr_res%status = ptr_res%status + NF_STATUS_NOT_CONVERGED

100 continue

    call result_update (ptr_res, x, fx, nit=iter, nfev=fcn%nfev)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)
end subroutine


subroutine __APPEND(root_halley_args,__PREC) (fcn, x, args, xtol, tol, &
        maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fcn_der2_args,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    real (PREC), intent(inout), dimension(:) :: args
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    call root_halley_impl (x, xtol, tol, maxiter, res, fcn_args=fcn, args=args)
end subroutine


subroutine __APPEND(root_halley,__PREC) (fcn, x, xtol, tol, maxiter, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fcn_der2,__PREC)) :: fcn
    real (PREC), intent(inout) :: x
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    call root_halley_impl (x, xtol, tol, maxiter, res, fcn=fcn)
end subroutine


subroutine __APPEND(root_halley_impl,__PREC) (x, xtol, tol, &
        maxiter, res, fcn, fcn_args, args)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(inout) :: x
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res
    procedure (__APPEND(fcn_der2,__PREC)), optional :: fcn
    procedure (__APPEND(fcn_der2_args,__PREC)), optional :: fcn_args
    real (PREC), intent(inout), dimension(:), optional :: args

    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res
    real (PREC) :: lxtol, ltol
    integer :: lmaxiter

    real (PREC) :: fx, fpx, fppx, x0, discr
    integer :: iter

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call check_inputs (xtol, tol, maxiter, res=ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    call set_defaults (xtol, tol, maxiter, lxtol, ltol, lmaxiter)

    x0 = x
    fppx = 0.0_PREC

    do iter = 1, lmaxiter
        if (present(fcn_args)) then
            call fcn_args (x0, args, fx, fpx, fppx)
        else
            call fcn (x0, fx, fpx, fppx)
        end if

        if (abs(fx) < ltol) then
            ptr_res%msg = "Convergence achieved; abs(f(x)) < ltol"
            ptr_res%status = NF_STATUS_OK
            goto 100
        end if

        if (fpx == 0.0_PREC) then
            ptr_res%msg = "Derivative evaluted to 0"
            ptr_res%status = NF_STATUS_OK
            goto 100
        end if

        ! Parabolic Halley's method
        discr = fpx ** 2.0_PREC - 2.0_PREC * fx * fppx
        if (discr < 0.0_PREC) then
            x = x0 - fpx / fppx
        else
            x = x0 - 2.0_PREC*fx / (fpx + signum(fpx) * sqrt(discr))
        end if

        ! Exit if tolerance level achieved
        if (abs(x - x0) < lxtol) then
            ptr_res%status = NF_STATUS_OK
            ptr_res%msg = "Convergence achieved: abs(x(n)-x(n-1)) < lxtol"
            goto 100
        end if

        ! otherwise update and go to next iteration
        x0 = x
    end do

   ptr_res%msg = "Max. number of iterations exceeded"
   ptr_res%status = NF_STATUS_MAX_ITER
   ptr_res%status = ptr_res%status + NF_STATUS_NOT_CONVERGED

100 continue

    call result_update (ptr_res, x, fx, nit=iter, nfev=iter)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

end subroutine
