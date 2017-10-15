program rosenbrock_slsqp

    use, intrinsic :: iso_fortran_env
    use numfort_optimize, optim_result => optim_result_real64

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()
    call example2 ()

contains

subroutine example1 ()

    real (PREC), parameter :: xlb(2) = [-1.0, -1.0], xub(2) = [1.0, 1.0]
    real (PREC), dimension(2) :: x0
    type (optim_result) :: res
    real (PREC), parameter :: tol = 1.0d-8

    integer :: m

    m = 1

    call minimize_slsqp (fobj, x0, xlb, xub, m, f_ieqcons=fconstr, tol=tol, &
        res=res)

    print '(tr1, "Min. at x=", 2(f0.8,:,", "), "; fx=", f0.8, "; iter=", i0, "; nfev=", i0)', &
        res%x, res%fx, res%nit, res%nfev

end subroutine


subroutine example2 ()
    real (PREC), parameter :: xlb(2) = [-1.0, -1.0], xub(2) = [1.0, 1.0]
    real (PREC), dimension(2) :: x0
    type (optim_result) :: res
    real (PREC), parameter :: tol = 1.0d-8

    integer :: m

    m = 1

    call minimize_slsqp (fobj, x0, xlb, xub, m, f_eqcons=fconstr2, tol=tol, &
        res=res)

    print '(tr1, "Min. at x=", 2(f0.8,:,", "), "; fx=", f0.8, "; iter=", i0, "; nfev=", i0)', &
        res%x, res%fx, res%nit, res%nfev


end subroutine

subroutine fobj (x, fx, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), optional, contiguous :: fpx

    if (present(fx)) then
        ! Compute objective
        fx = 1.0d2 * (x(2)-x(1)**2.0d0)**2.0d0 + (1.0d0-x(1))**2.0d0
    end if

    if (present(fpx)) then
        fpx(1) = - 4.0d2 * (x(2)-x(1)**2.0d0) * x(1) - 2.0d0*(1.0d0-x(1))
        fpx(2) = 2.0d2 * (x(2)-x(1)**2.0d0)
    end if
end subroutine


subroutine fconstr (x, fx, fpx)
    !*  Function evaluating inequality constraints
    !   Constraints needs to be formulated such that C(x) >= 0
    real (real64), intent(in), dimension(:), contiguous :: x
    real (real64), intent(out), dimension(:), contiguous, optional :: fx
    real (real64), intent(out), dimension(:,:), contiguous, optional :: fpx

    if (present(fx)) then
        fx = -x(1)**2.0d0 - x(2) ** 2.0d0 + 1.0d0
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
        fx = -x(1)**2.0d0 - x(2) ** 3.0d0 + 1.0d0
    end if

    ! Jacobian of constraint function
    if (present(fpx)) then
        fpx(1,1) = -2.0d0 * x(1)
        fpx(1,2) = -3.0d0 * x(2) ** 2.0d0
    end if
end subroutine


end program
