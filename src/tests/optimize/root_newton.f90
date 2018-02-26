

program test_optimize_root_newton

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_testing
    use numfort_optimize, only: root_newton, root_halley, &
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
    real (PREC), intent(out) :: fx, fpx
    call fcn_linear (x, fx)
    call jac_linear (x, fpx)
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
    real (PREC), intent(out) :: fx, fpx
    call fcn_linear_args (x, args, fx)
    call jac_linear_args (x, args, fpx)
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
    real (PREC), intent(out) :: fx, fpx
    call fcn_quad (x, fx)
    call jac_quad (x, fpx)
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
    real (PREC), intent(out) :: fx, fpx
    call fcn_quad_args (x, args, fx)
    call jac_quad_args (x, args, fpx)
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



subroutine fcn_jac_cycles (x, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx, fpx

    fx = x**3.0 - 2.0 * x + 2
    fpx = 3.0 * x**2.0 - 2.0
end subroutine

end program
