program newton_demo
    !*  Demo program for the 1-dimensional root finders using Newton's and
    !   Halley's methods

    use, intrinsic :: iso_fortran_env
    use numfort_optimize, optim_result => optim_result_real64

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()

contains


subroutine example1 ()

    real (PREC) :: x, x_true
    type (optim_result) :: res
    real (PREC), parameter :: tol = 1d-8, xtol = 1d-8

    ! Newton's method
    x = -10.0
    x_true = (-3.0 - sqrt(5.0d0)) / 2.0
    call root_newton (func1, x, xtol=xtol, tol=tol, res=res)
    call print_report (res, x_true)

    ! Halley's method using second derivatives
    x = -10.0
    call root_halley (func2, x, xtoL=xtol, tol=tol, res=res)
    call print_report (res, x_true)

end subroutine

subroutine func1 (x, fx, fpx, args)
    !   Wrapper for func2 that returns only f(x) and f'(x), but not the
    !   second derivative.
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx, fpx
    real (PREC), intent(in out), dimension(:), optional :: args

    real (PREC) :: fppx

    call func2 (x, fx, fpx, fppx, args)
end subroutine

subroutine func2 (x, fx, fpx, fppx, args)
    ! Objective function: 3rd-degree polynomial with roots at (-3-sqrt(5))/2,
    ! (-3+sqrt(5))/2 and 2
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx, fpx, fppx
    real (PREC), intent(in out), dimension(:), optional :: args

    ! use 3rd-degree polynomial with roots at (-3-sqrt(5))/2,
    ! (-3+sqrt(5))/2 and 2
    fx = (x+3) * (x-1)**2 - 5
    fpx = 3 * x ** 2 + 2 * x
    fppx = 6 * x + 2
end subroutine

subroutine print_report (res, exact_root)
    type (optim_result), intent(in) :: res
    real (PREC), intent(in) :: exact_root
    integer, save :: ii = 1

    print "('#', t3, 'Example ', i0)", ii
    print "(t3, 'Exit status: ', a)", char(res%status)
    if (len_trim(res%msg) > 0) then
        print "(t3, 'Message: ', a)", res%msg
    end if
    if (res%success) then
        print "(t3, 'Root located at: ', g23.15e2, '; true root: ', g23.15e2)", &
            res%x, exact_root
        print "(t3, 'Function value at root: ', g23.15e2)", res%fx
        print "(t3, 'Number of iterations: ', i0)", res%nit
        print "(t3, 'Number of function evaluations: ', i0)", res%nfev
    endif

    ii = ii + 1

end subroutine

end
