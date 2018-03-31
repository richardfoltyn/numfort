
program test_optimize_simplex
    !*  Unit tests for Simplex wrapper.

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env

    use fcore_common, FC_status_t => status_t
    use fcore_testing
    use numfort_arrays
    use numfort_common_testing
    use numfort_optimize, workspace => workspace_real64, &
        optim_result => optim_result_real64

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("Simplex wrapper unit tests")

    call test_quad (tests)

    call tests%print ()
end subroutine


subroutine test_quad (tests)
    class (test_suite), intent(inout) :: tests

    class (test_case), pointer :: tc
    type (optim_result) :: res
    integer, parameter :: NDIM = 4
    real (PREC), dimension(NDIM) :: x, x_desired

    tc => tests%add_test ("Quadratic objective function")

    x_desired = 0.0

    ! Initial starting point
    call arange(x)

    call minimize_simplex (fcn_quad, x, tol=1.0d-5, res=res)

    call tc%assert_true (res%status == NF_STATUS_OK, "Exit status")
    call tc%assert_true (all(res%x == x), "X == RES%X")
    call tc%assert_true (all_close (res%x, x_desired, atol=1.0d-3), &
        "Correct minimizer X")
    call tc%assert_true (abs(res%fx(1)) < 1.0d-3, &
        "Correct minimal function value f(X)")

end subroutine


subroutine fcn_quad (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx

    fx = sum(x**2.0_PREC)
end subroutine


end program