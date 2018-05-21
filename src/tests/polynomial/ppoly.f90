

program test_ppoly
    !*  Unit tests for function approximation with piecewise polynomials

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_common_testing
    use numfort_polynomial
    use numfort_interpolate

    use fcore_testing, only: test_suite, test_case
    use fcore_strings

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ('Piecewise polynomial unit tests')

    call test_ppolyfit_bernstein (tests)

    call tests%print ()

end subroutine



subroutine test_ppolyfit_bernstein (tests)
    !*  Unit tests for fitting and evaluating piecewise polynomials wrt. the
    !   Bernstein basis.
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(:), allocatable :: knots, coefs, x, p, xp, y1
    real (PREC), dimension(:,:), allocatable :: y
    integer :: k, ncoefs, nknots, nx, n
    real (PREC) :: xlb, xub
    type (ppoly_bernstein) :: pp
    type (status_t) :: status

    tc => tests%add_test ('PPOLYFIT tests using Bernstein basis')

    n = 51
    xlb = 1.0d-8
    xub = 100
    k = 3
    allocate (x(n), y(2, n))

    ! Compute function values and first derivatives
    call linspace (x, xlb, xub)
    call fcn2 (x, y(1,:), y(2,:))

    ncoefs = ppoly_get_ncoefs (pp, n, k)
    nknots = ppoly_get_nknots (pp, n, k)
    allocate (knots(nknots), coefs(ncoefs))

    ! Fit piecewise polynomial to 0th and 1st derivative data
    call ppolyfit (pp, x, y, k, knots, coefs, status=status)

    allocate (p(n))
    call ppolyval (pp, knots, coefs, x, p, status=status)
    call tc%assert_true (status == NF_STATUS_OK &
        .and. all_close (p, y(1,:), atol=1.0e-8_PREC), &
        'Fit / evaluate cubic piecewise polynomial at original knots')
    deallocate (p)

    ! Allocate at different points than knots
    nx = 107
    allocate (p(nx), xp(nx), y1(nx))
    call linspace (xp, xlb, xub)
    ! True function values
    call fcn2 (xp, y1)
    call ppolyval (pp, knots, coefs, xp, p, status=status)
    call tc%assert_true (status == NF_STATUS_OK &
        .and. all_close (p, y1, atol=1.0e-8_PREC), &
        'Evaluate cubic piecewise polynomial at other points')
    deallocate (p)


end subroutine


elemental subroutine fcn1 (x, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), optional :: fpx

    if (present(fx)) fx = log(x)
    if (present(fpx)) fpx = 1.0_PREC / x
end subroutine


elemental subroutine fcn2 (x, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), optional :: fpx

    real (PREC) :: lfx, lfpx
    real (PREC), parameter :: coefs(*) = [-123.34d0, 4.56d0, -7.23d0, 0.234d0]
    real (PREC) :: coefs_d(3)

    call polyval (coefs, x, lfx)
    call polyder (coefs, coefs_d, m=1)
    call polyval (coefs_d, x, lfpx)

    if (present(fx)) fx = lfx
    if (present(fpx)) fpx = lfpx
end subroutine



end