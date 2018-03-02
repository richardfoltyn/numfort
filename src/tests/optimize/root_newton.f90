

program test_optimize_root_newton

    use, intrinsic :: iso_fortran_env
    
    use numfort_core, only : PI_real64, signum
    use numfort_common
    use numfort_common_testing
    use numfort_optimize, only: root_newton, root_halley, root_newton_bisect, &
        optim_result => optim_result_real64

    use fcore_testing, only: test_suite, test_case
    use fcore_strings

    implicit none

    integer, parameter :: PREC = real64

    real (PREC), parameter :: LINEAR_CONST = -2.0
    real (PREC), parameter :: LINEAR_SLOPE = 2.0

    real (PREC), parameter :: QUAD_ROOT1 = -1.0
    real (PREC), parameter :: QUAD_ROOT2 = 1.0

    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("ROOT_NEWTON unit tests")

    call test_linear (tests)
    call test_quad (tests)
    call test_cycles (tests)
    call test_halley (tests)
    call test_newton_bisect (tests)

    call tests%print ()
end subroutine


subroutine test_linear (tests)
    type (test_suite), intent(in out) :: tests

    class (test_case), pointer :: tc

    real (PREC) :: x, fx
    real (PREC), dimension(2) :: args
    logical, parameter :: NDIFF = .true.
    type (optim_result) :: res
    real (PREC), parameter :: atol = 1.0d-12, rtol = 0.0
    real (PREC), parameter :: root_exact = 1.0, zero = 0.0
        !   True root and function value at root
    real (PREC), parameter :: x0 = 1.234

    tc => tests%add_test ("ROOT_NEWTON: Linear function")

    args(1) = LINEAR_CONST
    args(2) = LINEAR_SLOPE

    x = x0
    call root_newton (fcn_linear, x, NDIFF, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Num. derivative, no ARGS")

    x = x0
    call root_newton (fcn_linear, jac_linear, x, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Separate FCN and JAC, no ARGS")

    x = x0
    call root_newton (fcn_jac_linear, x, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Joint FCN and JAC, no ARGS")

    x = x0
    call root_newton (fcn_linear_args, x, args, NDIFF, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Num. derivate + ARGS")

    x = x0
    call root_newton (fcn_linear_args, jac_linear_args, x, args, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Separate FCN and JAC + ARGS")

    x = x0
    call root_newton (fcn_jac_linear_args, x, args, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Joint FCN and JAC + ARGS")


end subroutine



subroutine test_quad (tests)
    type (test_suite), intent(in out) :: tests

    class (test_case), pointer :: tc

    real (PREC) :: x, fx
    real (PREC), dimension(2) :: args
    logical, parameter :: NDIFF = .true.
    type (optim_result) :: res
    real (PREC), parameter :: atol = 1.0d-12, rtol = 0.0
    real (PREC), parameter :: root_exact = -1.0, zero = 0.0
        !   True root and function value at root
    real (PREC), parameter :: x0 = -2.234d0

    tc => tests%add_test ("ROOT_NEWTON: Quadratic function")

    args(1) = QUAD_ROOT1
    args(2) = QUAD_ROOT2

    x = x0
    call root_newton (fcn_quad, x, NDIFF, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Num. derivative, no ARGS")

    x = x0
    call root_newton (fcn_quad, jac_quad, x, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Separate FCN and JAC, no ARGS")

    x = x0
    call root_newton (fcn_jac_quad, x, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Joint FCN and JAC, no ARGS")

    x = x0
    call root_newton (fcn_quad_args, x, args, NDIFF, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Num. derivate + ARGS")

    x = x0
    call root_newton (fcn_quad_args, jac_quad_args, x, args, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Separate FCN and JAC + ARGS")

    x = x0
    call root_newton (fcn_jac_quad_args, x, args, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Joint FCN and JAC + ARGS")

end subroutine



subroutine test_newton_bisect (tests)
    type (test_suite), intent(in out) :: tests

    class (test_case), pointer :: tc

    real (PREC) :: x, fx
    type (optim_result) :: res
    real (PREC), dimension(1) :: args
    real (PREC), parameter :: atol = 1.0d-7, rtol = 0.0
    real (PREC), parameter :: zero = 0.0
    real (PREC) :: root_exact
        !   True root and function value at root
    real (PREC) :: x0

    tc => tests%add_test ("ROOT_NEWTON_BISECT unit tests")

    ! Check that error is raised the invalid bracket is passed
    x = 1.0
    call root_newton_bisect (fcn_jac_quad, x, a=1.0_PREC, b=1.0_PREC, res=res)
    call tc%assert_true (res%status == NF_STATUS_INVALID_ARG, &
        "Calling with invalid bracket a > b")

    call root_newton_bisect (fcn_jac_quad, x, a=2.0_PREC, b=3.0_PREC, res=res)
    call tc%assert_true (res%status == NF_STATUS_INVALID_ARG, &
        "Calling with invalid bracket x < a < b")

    call root_newton_bisect (fcn_jac_quad, x, a=-1.0_PREC, b=0.0_PREC, res=res)
    call tc%assert_true (res%status == NF_STATUS_INVALID_ARG, &
        "Calling with invalid bracket a < b < x")

    ! Test withr regular quadratic function which should work exactly as
    ! regular ROOT_NEWTON and converge in the same number of iterations.
    x0 = -2.234d0
    x = x0
    call root_newton_bisect (fcn_jac_quad, x, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, -1.0_PREC, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol) .and. res%success, &
        "Quadratic function, " // str(x0, 'f0.6'))

    x0 = 2.50d0
    root_exact = 0.0
    ! Test with regular Newton root finder
    x = x0
    call root_newton (fcn_jac_logistic, x, res=res)

    x = x0
    call root_newton_bisect (fcn_jac_logistic, x, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol) .and. res%success, &
        "Modified logistic CDF, " // str(x0, 'f0.6'))

    ! Test with regular Newton root-finder. This will fail for this
    ! particular starting point.
    x0 =  PI_real64 / 2.0 * 0.99d0
    root_exact = -0.739086_PREC
    x = x0
    call root_newton (fcn_jac_cosxx, x, res=res)
    
    x = x0
    call root_newton_bisect (fcn_jac_cosxx, x, res=res)
    fx = res%fx(1)
    ! True root according to Wolfram Alpha        
    call tc%assert_true (all_close (x, root_exact, atol=1.0e-5_PREC, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol) .and. res%success, &
        "cos(x)+x, " // str(x0, 'f0.6'))

    ! Function f(x) = |x|^alpha for alpha in (0,1/2)
    ! Standard ROOT_NEWTON will diverge
    x0 = -20.0d0
    root_exact = - 1.0 / (25 * sqrt(5.0_PREC))
    ! Pass alpha in ARGS array
    args(1) = 0.4d0

    x = x0
    call root_newton (fcn_jac_abs_sqrt, x, args, res=res)

    x = x0
    call root_newton_bisect (fcn_jac_abs_sqrt, x, args, a=-21.0d0, b=-0.01d0, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol) .and. res%success, &
        "|x|^alpha, " // str(x0, 'f0.6'))

end subroutine



subroutine test_halley (tests)
    type (test_suite), intent(in out) :: tests

    class (test_case), pointer :: tc

    real (PREC) :: x, fx
    real (PREC), dimension(2) :: args
    type (optim_result) :: res
    real (PREC), parameter :: atol = 1.0d-12, rtol = 0.0
    real (PREC), parameter :: root_exact = -1.0, zero = 0.0
        !   True root and function value at root
    real (PREC), parameter :: x0 = -2.234d0

    tc => tests%add_test ("ROOT_HALLEY: Quadratic function")

    args(1) = QUAD_ROOT1
    args(2) = QUAD_ROOT2

    x = x0
    call root_halley (fcn_halley, x, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Quadratic function, no ARGS")

    x = x0
    call root_halley (fcn_halley_args, x, args, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol), &
        "Quadratic function, ARGS")

end subroutine


subroutine test_cycles (tests)
    type (test_suite), intent(in out) :: tests

    class (test_case), pointer :: tc

    real (PREC) :: x, fx
    type (optim_result) :: res
    real (PREC), parameter :: atol = 1.0d-7, rtol = 0.0
    real (PREC), parameter :: root_exact = -1.76929235d0, zero = 0.0
        !   True root and function value at root
    real (PREC), parameter :: x0 = 0.99d0

    tc => tests%add_test ("ROOT_NEWTON: period-2 cycles")

    ! Check objective function from wikipedia that is known to produce
    ! cycles if x0 is close to 1.0 or 0.0.

    x = x0
    call root_newton (fcn_jac_cycles, x, xtol2=1.0e-8_PREC, maxiter=100, res=res)
    fx = res%fx(1)
    call tc%assert_true (all_close (x, root_exact, atol=atol, rtol=rtol) .and. &
        all_close (fx, zero, atol=atol, rtol=rtol) .and. res%success, &
        "f(x) = x**3.0 - 2*x + 2; x(0)=" // str(x0, 'f0.6'))

end subroutine



subroutine fcn_linear (x, fx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    fx = LINEAR_CONST + LINEAR_SLOPE * x
end subroutine

subroutine jac_linear (x, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fpx
    fpx = LINEAR_SLOPE
end subroutine

subroutine fcn_jac_linear (x, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx, fpx
    
    if (present(fx)) then
        call fcn_linear (x, fx)
    end if
    
    if (present(fpx)) then
        call jac_linear (x, fpx)
    end if
end subroutine

subroutine fcn_linear_args (x, args, fx)
    real (PREC), intent(in) :: x
    real (PREC), intent(in out), dimension(:) :: args
    real (PREC), intent(out) :: fx
    fx = args(1) + args(2) * x
end subroutine

subroutine jac_linear_args (x, args, fx)
    real (PREC), intent(in) :: x
    real (PREC), intent(in out), dimension(:) :: args
    real (PREC), intent(out) :: fx
    fx = args(2)
end subroutine

subroutine fcn_jac_linear_args (x, args, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(in out), dimension(:) :: args
    real (PREC), intent(out), optional :: fx, fpx
    
    if (present(fx)) then
        call fcn_linear_args (x, args, fx)
    end if
    
    if (present(fx)) then
        call jac_linear_args (x, args, fpx)
    end if
end subroutine



subroutine fcn_quad (x, fx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    fx = (x-QUAD_ROOT1) * (x-QUAD_ROOT2)
end subroutine

subroutine jac_quad (x, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fpx
    fpx = 2.0 * x - (QUAD_ROOT1 + QUAD_ROOT2)
end subroutine

subroutine fcn_jac_quad (x, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx, fpx
    
    if (present(fx)) then
        call fcn_quad (x, fx)
    end if
    
    if (present(fpx)) then
        call jac_quad (x, fpx)
    end if
end subroutine

subroutine fcn_quad_args (x, args, fx)
    real (PREC), intent(in) :: x
    real (PREC), intent(in out), dimension(:) :: args
    real (PREC), intent(out) :: fx
    fx = (x-args(1)) * (x-args(2))
end subroutine

subroutine jac_quad_args (x, args, fx)
    real (PREC), intent(in) :: x
    real (PREC), intent(in out), dimension(:) :: args
    real (PREC), intent(out) :: fx
    fx = 2.0 * x - (sum(args))
end subroutine

subroutine fcn_jac_quad_args (x, args, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(in out), dimension(:) :: args
    real (PREC), intent(out), optional :: fx, fpx
    
    if (present(fx)) then
        call fcn_quad_args (x, args, fx)
    end if
    
    if (present(fpx)) then
        call jac_quad_args (x, args, fpx)
    end if
end subroutine



subroutine fcn_halley (x, fx, fpx, fppx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx, fpx, fppx
    fx = (x-QUAD_ROOT1) * (x-QUAD_ROOT2)
    fpx = 2.0 * x - (QUAD_ROOT1 + QUAD_ROOT2)
    fppx = 2.0
end subroutine

subroutine fcn_halley_args (x, args, fx, fpx, fppx)
    real (PREC), intent(in) :: x
    real (PREC), intent(in out), dimension(:) :: args
    real (PREC), intent(out) :: fx, fpx, fppx
    fx = (x-args(1)) * (x-args(2))
    fpx = 2.0 * x - (args(1) + args(2))
    fppx = 2.0
end subroutine



subroutine fcn_jac_cycles (x, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx, fpx

    if (present(fx)) then
        fx = x**3.0 - 2.0 * x + 2
    end if
    
    if (present(fpx)) then
        fpx = 3.0 * x**2.0 - 2.0
    end if
end subroutine



subroutine fcn_jac_logistic (x, fx, fpx)
    !*  Modified logistic CDF
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx, fpx

    real (PREC) :: ex1

    ex1 = exp(-x) + 1.0_PREC

    if (present(fx)) then
        fx = 1.0_PREC/ex1 - 0.5_PREC
    end if
    
    if (present(fpx)) then
        fpx = exp(-x) / ex1**2.0_PREC
    end if
end subroutine


subroutine fcn_jac_cosxx (x, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx, fpx

    if (present(fx)) then
        fx = cos(x) + x
    end if
    
    if (present(fpx)) then
        fpx = -sin(x) + 1.0_PREC
    end if
end subroutine


subroutine fcn_jac_abs_sqrt (x, args, fx, fpx)
    !*  Function f(x) = |x|^alpha.
    !   Alpha should be (0,1/2) for overshooting and divergence with regular
    !   Newton method.
    real (PREC), intent(in) :: x
    real (PREC), intent(in out), dimension(:) :: args
    real (PREC), intent(out), optional :: fx, fpx

    if (present(fx)) then
        fx = abs(x) ** args(1) - 0.2
    end if
    
    ! Ignore non-differentiable point x = 0
    if (present(fpx)) then
        fpx = signum (x) * args(1) * abs(x) ** (args(1) - 1.0_PREC)
    end if

end subroutine

end program
