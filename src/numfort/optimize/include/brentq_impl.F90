
subroutine __APPEND(root_brentq,__PREC) (f, a, b, xtol, rtol, maxiter, x0, res)

    integer, parameter :: PREC = __PREC

    procedure (__APPEND(fcn,__PREC)) :: f
    real (PREC), intent(in) :: a, b
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: rtol
    real (PREC), intent(out) :: x0
    integer, intent(in), optional :: maxiter
    type (__APPEND(optim_result,__PREC)), intent(in out), optional :: res

    real (PREC) :: lxtol, lrtol, fx0
    ! maximum iterations, number of function evaluations
    integer :: lmaxiter, nfev, iter
    type (status_t) :: status

    lmaxiter = 100
    lrtol = sqrt(epsilon(1.0_PREC))
    lxtol = min(1.0e4_PREC * epsilon(1.0_PREC), 1.0e-5_PREC)

    if (present(maxiter)) lmaxiter = maxiter
    if (present(xtol)) lxtol = xtol
    if (present(rtol)) lrtol = rtol

    call root_brentq_impl (f, a, b, lxtol, lrtol, lmaxiter, x0, fx0, iter, nfev, status)

    if (present(res)) then
        if (NF_STATUS_INVALID_ARG .in. status) then
            call result_update (res, status=status, msg="Not a bracketing interval")
        else
            call result_update (res, x=x0, fx=fx0, status=status, nit=iter, nfev=nfev)
        end if
    end if

end subroutine

subroutine __APPEND(root_brentq_impl,__PREC) (f, a, b, xtol, rtol, maxiter, &
        x0, fx0, iter, nfev, status)
    !*  ROOT_BRENTQ_IMPL is a F90 port of Brent's method as implemented
    !   in Scipy's brentq()

    integer, parameter :: PREC = __PREC

    procedure (__APPEND(fcn,__PREC)) :: f
    real (PREC) :: a, b, xtol, rtol, x0, fx0
    integer :: maxiter, nfev, iter
    type (status_t), intent(out) :: status

    intent (in) :: a, b, xtol, rtol, maxiter
    intent (out) :: iter, nfev, x0, fx0

    real (PREC) :: xpre, xcur, xblk, fcur, fpre, fblk, spre, scur, sbis
    real (PREC) :: delta, stry, dpre, dblk

    integer :: i

    xpre = a
    xcur = b
    xblk = 0.0_PREC
    fblk = 0.0_PREC
    spre = 0.0_PREC
    scur = 0.0_PREC

    call f (xpre, fpre)
    call f (xcur, fcur)
    nfev = 2
    iter = 0

    if (fpre * fcur > 0.0_PREC) then
        status = NF_STATUS_INVALID_ARG
        return
    else if (fpre == 0.0_PREC) then
        status = NF_STATUS_OK
        x0 = xpre
        fx0 = fpre
        return
    else if (fcur == 0.0_PREC) then
        status = NF_STATUS_OK
        x0 = xcur
        fx0 = fcur
        return
    end if

    do i = 1, maxiter
        iter = iter + 1

        if (fpre * fcur < 0) then
            xblk = xpre
            fblk = fpre
            spre = xcur - xpre
            scur = spre
        end if

        ! make sure that xcur is the better candidate root, ie f(xcur) is closer
        ! to zero; if needed, swap xcur / xblk
        if (abs(fblk) < abs(fcur)) then
            xpre = xcur
            xcur = xblk
            xblk = xpre

            fpre = fcur
            fcur = fblk
            fblk = fpre
        end if

        delta = xtol + rtol * abs(xcur) / 2.0_PREC
        sbis = (xblk - xcur) / 2.0_PREC

        ! check for convergence
        if (fcur == 0.0_PREC .or. abs(sbis) < delta) then
            x0 = xcur
            fx0 = fcur
            status = NF_STATUS_OK
            return
        end if

        ! compute step size which determines updated xcur
        if (abs(spre) > delta .and. abs(fcur) < abs(fpre)) then
            if (xpre == xblk) then
                ! secant method
                stry = -fcur * (xcur - xpre) / (fcur - fpre)
            else
                ! extrapolate
                dpre = (fpre - fcur) / (xpre - xcur)
                dblk = (fblk - fcur) / (xblk - xcur)
                stry = -fcur * (fblk*dblk - fpre*dpre) / (dblk*dpre*(fblk - fpre))
            end if

            if (2 * abs(stry) < min(abs(spre), 3.0*abs(sbis) - delta)) then
                ! good short step
                spre = scur
                scur = stry
            else
                ! bisect
                spre = sbis
                scur = sbis
            end if
        else
            ! bisect
            spre = sbis
            scur = sbis
        end if

        xpre = xcur
        fpre = fcur

        if (abs(scur) > delta) then
            xcur = xcur + scur
        else
            if (sbis > 0) then
                xcur = xcur + delta
            else
                xcur = xcur - delta
            end if
        end if

        call f (xcur, fcur)
        nfev = nfev + 1
    end do

    ! algorithm did not converge
    status = NF_STATUS_NOT_CONVERGED
    ! store best guesses so far
    x0 = xcur
    fx0 = fcur

end subroutine
