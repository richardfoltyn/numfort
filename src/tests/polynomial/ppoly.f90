

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

    call test_ppolyfit_bernstein_input (tests)
    call test_ppolyfit_bernstein (tests)

    call tests%print ()

end subroutine



subroutine test_ppolyfit_bernstein (tests)
    !*  Unit tests for fitting and evaluating piecewise polynomials wrt. the
    !   Bernstein basis.
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(:), allocatable :: knots, coefs, coefs1, x, p, xp, y1
    real (PREC), dimension(:,:), allocatable :: y
    integer :: k, ncoefs, nknots, nx, n
    real (PREC) :: xlb, xub, left, right
    type (ppoly_bernstein) :: pp, pp1
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
    call tc%assert_true (status == NF_STATUS_OK, &
        'Cubic: fit piecewise polynomial')

    allocate (p(n))
    call ppolyval (pp, knots, coefs, x, p, status=status)
    call tc%assert_true (status == NF_STATUS_OK &
        .and. all_close (p, y(1,:), atol=1.0e-8_PREC), &
        'Cubic: evaluate at original knots')
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
        'Cubic: Evaluate at points other than knots')
    deallocate (p, xp, y1)

    ! Extrapolation
    nx = 17
    allocate (p(nx), xp(nx), y1(nx))
    call linspace(xp, xlb-1.0, xub+1.0)
    ! Some additional non-interior points
    xp(5) = xlb - 10.0_PREC
    xp(15) = xub + 10.0_PREC
    ! Correct function values
    call fcn2 (xp, y1)
    call ppolyval (pp, knots, coefs, xp, p, ext=NF_INTERP_EVAL_EXTRAPOLATE, &
        status=status)
    call tc%assert_true (status == NF_STATUS_OK &
        .and. all_close (p, y1, atol=1.0e-8_PREC), &
        'Cubic: Evaluate at non-interior points')

    ! Replace non-interior points with come constant values
    left = -10000.0
    right = 10000.0
    where (xp < xlb)
        y1 = left
    else where (xp > xub)
        y1 = right
    end where
    call ppolyval (pp, knots, coefs, xp, p, ext=NF_INTERP_EVAL_CONST, &
        left=left, right=right, status=status)
    call tc%assert_true (status == NF_STATUS_OK &
        .and. all_close (p, y1, atol=1.0e-8_PREC), &
        'Cubic: Evaluate, assign const to non-interior points')

    ! Raise error on encountering non-interior points
    call ppolyval (pp, knots, coefs, xp, p, ext=NF_INTERP_EVAL_ERROR, &
        status=status)
    call tc%assert_true (status == NF_STATUS_BOUNDS_ERROR, &
        'Cubic: Evaluate, raise error on non-interior points')
    deallocate (xp, y1, p)

    ! Test evaluating first derivative
    ncoefs = ppoly_get_ncoefs (pp, k=k-1)
    allocate (coefs1(ncoefs))
    call ppolyder (pp, knots, coefs, pp1, coefs1, m=1, status=status)

    nx = 13
    status = NF_STATUS_UNDEFINED
    allocate (xp(nx), y1(nx), p(nx))
    call linspace (xp, xlb, xub)
    ! True first derivative
    call fcn2 (xp, fpx=y1)
    call ppolyval (pp1, knots, coefs1, xp, p, status=status)
    call tc%assert_true (status == NF_STATUS_OK &
        .and. all_close (p, y1, atol=1.0e-8_PREC), &
        'Cubic: Evaluate first derivative')

    ! Extrapolate first derivative
    status = NF_STATUS_UNDEFINED
    call linspace (xp, xlb-10.0, xub+10.0)
    call fcn2 (xp, fpx=y1)
    ! Note: extrapolates by default
    call ppolyval (pp1, knots, coefs1, xp, p, status=status)
    call tc%assert_true (status == NF_STATUS_OK &
        .and. all_close (p, y1, atol=1.0e-8_PREC), &
        'Cubic: Evaluate first derivative at non-interior points')
    deallocate (xp, p, y1)
end subroutine


subroutine test_ppolyfit_bernstein_input (tests)
    !*  Unit tests for input checking for Bernstein basis fitting routines
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(:), allocatable :: knots, coefs
    real (PREC), dimension(3) :: x
    real (PREC), dimension(2,3) :: y
    integer :: n, k, ncoefs, nknots
    type (ppoly_bernstein) :: pp
    type (status_t) :: status

    tc => tests%add_test ('PPOLYFIT input checks for Bernstein basis')

    ! Incompatible data point arrays
    k = 3
    x = 0.0
    y = 0.0
    n = size(x)
    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))
    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x(1:2), y, k, knots, coefs, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Incompatible number of data points')

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x, y(1:1,:), k, knots, coefs, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Too few derivatives for requested polynomial degree')

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x, y(:0,:), -1, knots, coefs, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid polynomial degree')

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x, y, k, knots(1:nknots-1), coefs, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Incompatible KNOTS array size')

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x, y, k, knots, coefs(1:ncoefs-1), status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Incompatible COEFS array size')
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