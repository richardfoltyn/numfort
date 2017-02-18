module numfort_optimize_brent

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_optim_result_mod

    implicit none
    private

    integer, parameter :: PREC = real64

    interface
        function func (x) result(fx)
            import PREC
            real (PREC), intent(in) :: x
            real (PREC) :: fx
        end function
    end interface

    public :: brentq


contains

subroutine brentq (f, a, b, xtol, rtol, maxiter, x0, res)

    procedure (func) :: f
    real (PREC), intent(in) :: a, b, xtol, rtol
    real (PREC), intent(out) :: x0
    integer, intent(in) :: maxiter
    class (optim_result), intent(in out), optional :: res

    optional :: maxiter, xtol, rtol

    real (PREC) :: lxtol, lrtol, fx0, lx0(1)
    ! maximum iterations, number of function evaluations
    integer :: lmaxiter, nfev, iter, status

    lmaxiter = 100
    lrtol = 6.0 * epsilon(1.0_PREC)
    lxtol = 2d-12

    if (present(maxiter)) lmaxiter = maxiter
    if (present(xtol)) lxtol = xtol
    if (present(rtol)) lrtol = rtol

    call brentq_impl (f, a, b, lxtol, lrtol, lmaxiter, x0, fx0, iter, nfev, status)

    if (present(res)) then
        if (status == NF_STATUS_INVALID_ARG) then
            call res%update (status=status, msg="Not a bracketing interval")
        else
            lx0(1) = x0
            call res%update (x=lx0, fx=fx0, status=status, nit=iter, nfev=nfev)
        end if
    end if


end subroutine

! BRENTQ is a F90 port of Brent's method as implemented in Scipy's brentq()
subroutine brentq_impl (f, a, b, xtol, rtol, maxiter, x0, fx0, iter, nfev, status)
    procedure (func) :: f
    real (PREC) :: a, b, xtol, rtol, x0, fx0
    integer :: maxiter, nfev, iter, status

    intent (in) :: a, b, xtol, rtol, maxiter
    intent (out) :: iter, nfev, status, x0, fx0

    real (PREC) :: xpre, xcur, xblk, fcur, fpre, fblk, spre, scur, sbis
    real (PREC) :: delta, stry, dpre, dblk

    integer :: i

    xpre = a
    xcur = b
    xblk = 0.0d0
    fblk = 0.0d0
    spre = 0.0d0
    scur = 0.0d0

    fpre = f (xpre)
    fcur = f (xcur)
    nfev = 2
    iter = 0

    if (fpre * fcur > 0.0) then
        status = NF_STATUS_INVALID_ARG
        return
    else if (fpre == 0.0) then
        status = NF_STATUS_OK
        return
    else if (fcur == 0.0) then
        status = NF_STATUS_OK
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

        fcur = f (xcur)
        nfev = nfev + 1
    end do

    ! algorithm did not converge
    status = NF_STATUS_NOT_CONVERGED
    ! store best guesses so far
    x0 = xcur
    fx0 = fcur

end subroutine

end module
