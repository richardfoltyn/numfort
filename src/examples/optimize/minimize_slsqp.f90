program rosenbrock_slsqp

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env
    use numfort_optimize, optim_result => optim_result_real64
    use numfort_common_workspace, workspace => workspace_real64

    implicit none

    integer, parameter :: PREC = real64

    integer :: nfev, ncev

    call example1 ()
    call example2 ()
    call example3 ()

contains

subroutine example1 ()

    real (PREC), parameter :: xlb(2) = [-1.0, -1.0], xub(2) = [1.0, 1.0]
    real (PREC), dimension(2) :: x0, x
    type (optim_result) :: res, res_ng
    real (PREC), parameter :: tol = 1.0d-8
    type (workspace) :: work
    real (PREC) :: fx
    integer, parameter :: m = 1
    real (PREC), dimension(m) :: c


    nfev = 0
    ncev = 0
    x0 = xlb + 0.5*(xub-xlb)
    x = x0
    call minimize_slsqp (fobj, x, xlb, xub, m, f_ieqcons=fconstr, tol=tol, &
        res=res)
    call fobj (x, fx)
    call fconstr (x, c)
    print 100, 'SLSQP', fx, c

    nfev = 0
    ncev = 0
    x = x0
    call minimize_slsqp_ng (fobj, x, work, fconstr, m, tol=tol, lbounds=xlb, &
        ubounds=xub, res=res_ng)
    call fobj (x, fx)
    call fconstr (x, c)
    print 100, 'SLSQP-NG', fx, c

    ! Repeat with exact line search
    x0 = xlb + 0.5*(xub-xlb)
!    call minimize_slsqp (fobj, x0, xlb, xub, m, f_ieqcons=fconstr, tol=tol, &
!        res=res, linesearch=NF_LINESEARCH_EXACT)

100 format (a10, ': f(x): ', es20.12e2, '; c(x): ', *(es20.12e2, :, ', '))

end subroutine



subroutine example2 ()
    real (PREC), parameter :: xlb(2) = [-1.0, -1.0], xub(2) = [1.0, 1.0]
    real (PREC), dimension(2) :: x0, x
    type (optim_result) :: res, res_ng
    real (PREC), parameter :: tol = 1.0d-8
    type (workspace) :: work
    real (PREC) :: fx
    integer, parameter :: m = 1
    real (PREC), dimension(m) :: c

    integer :: meq

    x0 = sqrt(0.5)
    meq = m

    ! Pick initial guess that satisfies equality constraint
    nfev = 0
    ncev = 0
    x = x0
    call minimize_slsqp (fobj, x, xlb, xub, m, f_eqcons=fconstr2, tol=tol, &
        res=res)
    call fobj(x, fx)
    call fconstr2 (x, c)
    print 100, 'SLSQP', fx, c

    nfev = 0
    ncev = 0
    x = x0
    call minimize_slsqp_ng (fobj, x, work, fconstr2, m, meq, lbounds=xlb, &
        ubounds=xub, tol=tol, res=res_ng)
    call fobj (x, fx)
    call fconstr2 (x, c)
    print 100, 'SLSQP-NG', fx, c

    ! Repeat with exact line search
    x0 = sqrt(0.5)
!    call minimize_slsqp (fobj, x0, xlb, xub, m, f_eqcons=fconstr2, tol=tol, &
!        res=res, linesearch=NF_LINESEARCH_EXACT)

100 format (a10, ': f(x): ', es20.12e2, '; c(x): ', *(es20.12e2, :, ', '))

end subroutine



subroutine example3 ()
    !*  Same as example 1, but without upper bounds

    real (PREC), parameter :: xlb(2) = [-1.0, -1.0]
    real (PREC), dimension(2) :: x0, xub, x
    type (optim_result) :: res, res_ng
    real (PREC), parameter :: tol = 1.0d-8
    real (PREC) :: POS_INF
    type (workspace) :: work
    integer, parameter :: m = 1
    real (PREC), dimension(m) :: c
    real (PREC) :: fx

    POS_INF = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    xub = POS_INF

    x0 = xlb + 0.1

    nfev = 0
    ncev = 0
    x = x0
    call minimize_slsqp (fobj, x, xlb, xub, m, f_ieqcons=fconstr, tol=tol, &
        res=res)
    call fobj (x, fx)
    call fconstr (x, c)
    print 100, 'SLSQP', fx, c

    nfev = 0
    ncev = 0
    x = x0
    call minimize_slsqp_ng (fobj, x, work, fconstr, m, lbounds=xlb, &
        ubounds=xub, tol=tol, res=res_ng)
    call fobj (x, fx)
    call fconstr (x, c)
    print 100, 'SLSQP-NG', fx, c

    ! Repeat with exact line search
    x0 = xlb + 0.1
!    call minimize_slsqp (fobj, x0, xlb, xub, m, f_ieqcons=fconstr, tol=tol, &
!        res=res, linesearch=NF_LINESEARCH_EXACT)

100 format (a10, ': f(x): ', es20.12e2, '; c(x): ', *(es20.12e2, :, ', '))

end subroutine



subroutine fobj (x, fx, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), optional, contiguous :: fpx

    if (present(fx)) then
        nfev = nfev + 1
        ! Compute objective
        fx = 1.0d2 * (x(2)-x(1)**2.0d0)**2.0d0 + (1.0d0-x(1))**2.0d0
!        print '(i3, " f(x): ", es16.8e2)', nfev, fx
    end if

    if (present(fpx)) then
        fpx(1) = - 4.0d2 * (x(2)-x(1)**2.0d0) * x(1) - 2.0d0*(1.0d0-x(1))
        fpx(2) = 2.0d2 * (x(2)-x(1)**2.0d0)
    end if
end subroutine



subroutine fconstr (x, fx, fpx)
    !*  Function evaluating inequality constraints
    !   Constraints needs to be formulated such that C(x) >= 0
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx

    if (present(fx)) then
        ncev = ncev + 1
        fx = -x(1)**2.0d0 - x(2) ** 2.0d0 + 1.0d0
!        print '(i3, " c(x): ", es16.8e2)', ncev, fx
    end if

    ! Jacobian of constraint function
    if (present(fpx)) then
        fpx(1,:) = -2.0d0 * x
    end if
end subroutine


subroutine fconstr2 (x, fx, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx

    if (present(fx)) then
        ncev = ncev + 1
        fx = x(1)**2.0d0 + x(2) ** 2.0d0 - 1.0d0
!        print '(i3, " c(x): ", es16.8e2)', ncev, fx
    end if

    ! Jacobian of constraint function
    if (present(fpx)) then
        fpx(1,1) = 2.0d0 * x(1)
        fpx(1,2) = 2.0d0 * x(2)
    end if
end subroutine


end program
