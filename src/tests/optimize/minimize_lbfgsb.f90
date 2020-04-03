
program test_optimize_minimize_lbfgsb

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common_testing
    use numfort_optimize, workspace => workspace_real64, &
        optim_result => optim_result_real64

    use fcore_common, FC_status_t => status_t
    use fcore_testing

    implicit none

    integer, parameter :: PREC = real64

    real (PREC), parameter :: ROOTS(*) = [0.123d0, -12.345d0, 1.2345d0]

    call test_all ()

    contains



subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("L-BFGS-B unit test")

    call test_quad (tests)

    call tests%print ()

end subroutine


subroutine test_quad (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: NDIM = 3
    type (workspace) :: work
    type (optim_result) :: res
    real (PREC), dimension(NDIM) :: x0
    logical, parameter :: ndiff = .true.

    tc => tests%add_test ("Quadratic objective")

    ! unconstrained minimizer
    x0 = 0.0
    call minimize_lbfgsb (fobj, x0, ndiff=ndiff, work=work, res=res)

    call tc%assert_true (all_close (res%x, ROOTS, atol=1.0d-6), &
        "Using numerical differentiation")

end subroutine


subroutine fobj (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx

    fx = (x(1)-ROOTS(1))**2.0d0 + (x(2)-ROOTS(2)) ** 2.0d0 + (x(3)-ROOTS(3))**2.0

end subroutine


end
