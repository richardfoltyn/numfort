

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

    call test_ppolyval_power_cubic (tests)

    call tests%print ()

end subroutine



subroutine test_ppolyfit_bernstein (tests)
    !*  Unit tests for fitting and evaluating piecewise polynomials wrt. the
    !   Bernstein basis.
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(:), allocatable :: knots, coefs, coefs1, xx, xp, yp, yy
    real (PREC), dimension(:), allocatable :: fx, fpx
    real (PREC), dimension(:,:), allocatable :: ydat
    integer :: k, ncoefs, nknots, nx, n, i
    real (PREC) :: xlb, xub, left, right, ylb, yub
    type (ppoly_bernstein) :: pp, pp1
    logical :: status_ok, values_ok
    type (status_t) :: status

    tc => tests%add_test ('PPOLYFIT tests using Bernstein basis')

    n = 11
    xlb = 1.0d-8
    xub = 100
    k = 3
    allocate (xx(n), ydat(2, n))
    allocate (fx(n), fpx(n))

    ! Compute function values and first derivatives
    call linspace (xx, xlb, xub)
    call fcn2 (xx, fx, fpx)
    ydat(1,:) = fx
    ydat(2,:) = fpx

    ncoefs = ppoly_get_ncoefs (pp, n, k)
    nknots = ppoly_get_nknots (pp, n, k)
    allocate (knots(nknots), coefs(ncoefs))

    ! Fit piecewise polynomial to 0th and 1st derivative data
    call ppolyfit (pp, xx, ydat, k, knots, coefs, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        '1d API: Cubic: fit piecewise polynomial')

    allocate (yy(n))
    call ppolyval (pp, knots, coefs, xx, yy, status=status)
    values_ok = all_close (yy, fx, atol=1.0e-8_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: Cubic: evaluate at original knots')
    deallocate (yy)

    ! Evaluate using scalar interface
    allocate (yy(n))
    status_ok = .true.
    do i = 1, size(xx)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xx(i), yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, fx, atol=1.0e-8_PREC)
    call tc%assert_true (status_ok .and. values_ok, &
        '0d API: Cubic: evaluate at original knots')
    deallocate (yy)

    ! === Evaluate at different points than knots ===
    nx = 107
    allocate (yp(nx), xp(nx), yy(nx))
    call linspace (xp, xlb, xub)
    ! True function values
    call fcn2 (xp, yp)
    call ppolyval (pp, knots, coefs, xp, yy, status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: Cubic: Evaluate at points other than knots')
    deallocate (yy)

    ! Evaluate using scalar interface
    allocate (yy(nx), source=0.0_PREC)
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xp(i), yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status_ok .and. values_ok,  &
        '0d API: Cubic: Evaluate at points other than knots')
    deallocate (yp, xp, yy)

    ! === Extrapolation ===
    nx = 17
    allocate (yp(nx), xp(nx), yy(nx))
    call linspace(xp, xlb-1.0, xub+1.0)
    ! Some additional non-interior points
    xp(5) = xlb - 10.0_PREC
    xp(15) = xub + 10.0_PREC
    ! Correct function values
    call fcn2 (xp, yp)
    call ppolyval (pp, knots, coefs, xp, yy, ext=NF_INTERP_EVAL_EXTRAPOLATE, &
        status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '0d API: Cubic: Evaluate at non-interior points')
    deallocate (yy)

    ! Scalar interface
    allocate (yy(nx), source=0.0_PREC)
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xp(i), yy(i), &
            ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status_ok .and. values_ok, &
        '0d API: Cubic: Evaluate at non-interior points')

    ! === Replace non-interior points with some constant values ===
    left = -10000.0
    right = 10000.0
    where (xp < xlb)
        yp = left
    else where (xp > xub)
        yp = right
    end where
    call ppolyval (pp, knots, coefs, xp, yy, ext=NF_INTERP_EVAL_CONST, &
        left=left, right=right, status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '0d API: Cubic: Evaluate, assign const to non-interior points')
    deallocate (yy)

    ! Scalar interface
    allocate (yy(nx))
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xp(i), yy(i), ext=NF_INTERP_EVAL_CONST, &
            left=left, right=right, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status_ok .and. values_ok, &
        '0d API: Cubic: Evaluate, assign const to non-interior points')

    ! === Raise error on encountering non-interior points ===
    call ppolyval (pp, knots, coefs, xp, yy, ext=NF_INTERP_EVAL_ERROR, &
        status=status)
    call tc%assert_true (status == NF_STATUS_BOUNDS_ERROR, &
        '1d API: Cubic: Evaluate, raise error on non-interior points')

    ! Replace non-interior points with boundary values
    call fcn2 (xp, yp)
    call fcn2 (xlb, ylb)
    call fcn2 (xub, yub)
    where (xp < xlb)
        yp = ylb
    else where (xp > xub)
        yp = yub
    end where

    call ppolyval (pp, knots, coefs, xp, yy, ext=NF_INTERP_EVAL_BOUNDARY, &
        status=status)
    values_ok = all_close (yp, yy, atol=1.0e-8_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: Cubic: Evaluate, assign const to non-interior points')

    deallocate (xp, yy, yp)

    ! === Test evaluating first derivative ===
    ncoefs = ppoly_get_ncoefs (pp, k=k-1)
    allocate (coefs1(ncoefs))
    call ppolyder (pp, knots, coefs, pp1, coefs1, m=1, status=status)

    nx = 13
    status = NF_STATUS_UNDEFINED
    allocate (xp(nx), yy(nx), yp(nx))
    call linspace (xp, xlb, xub)
    ! True first derivative
    call fcn2 (xp, fpx=yp)
    call ppolyval (pp1, knots, coefs1, xp, yy, status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: Cubic: Evaluate first derivative')
    deallocate (yy)

    ! Scalar interface
    allocate (yy(nx), source=0.0_PREC)
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp1, knots, coefs1, xp(i), yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status_ok .and. values_ok, &
        '0d API: Cubic: Evaluate first derivative')

    ! === Extrapolate first derivative ===
    status = NF_STATUS_UNDEFINED
    call linspace (xp, xlb-10.0, xub+10.0)
    call fcn2 (xp, fpx=yp)
    ! Note: extrapolates by default
    call ppolyval (pp1, knots, coefs1, xp, yy, status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: Cubic: Evaluate first derivative at non-interior points')
    deallocate (yy)

    ! Scalar interface
    allocate (yy(nx), source=0.0_PREC)
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp1, knots, coefs1, xp(i), yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    call tc%assert_true (status_ok .and. values_ok, &
        '0d API: Cubic: Evaluate first derivaitve at non-interior points')

    deallocate (xp, yp, yy)

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



subroutine test_ppolyval_power_cubic (tests)
    !*  Unit tests for evaluating cubic polynomials and their derivatives
    !   using a power basis.
    class (test_suite) :: tests

    type (ppoly) :: pp, pp_d1
    type (ppoly_bernstein) :: pp_bern, pp_bern_d1

    real (PREC), dimension(:), allocatable :: xx, yy, yy_bern
    real (PREC), dimension(:,:), allocatable :: ydat
    real (PREC), dimension(:), allocatable :: knots, coefs_bern, coefs
    real (PREC), dimension(:), allocatable :: coefs_bern_d1, coefs_d1
    real (PREC) :: lb, ub
    integer :: n, deg, nknots, ncoefs, ncoefs_bern
    type (status_t) :: status
    logical :: is_ok

    class (test_case), pointer :: tc

    tc => tests%add_test ('Unit tests for eval. cubic polynomials (power basis)')

    ! Fit some polynomial using the Bernstein basis
    nknots = 20
    deg = 3
    allocate (xx(nknots), ydat(2,nknots))

    lb = 0.0
    ub = 10.0
    call linspace (xx, 0.0_PREC, 10.0_PREC)
    ydat(1,:) = sin(xx)
    ydat(2,:) = cos(xx)

    ncoefs = ppoly_get_ncoefs (pp_bern, n=nknots, k=deg)
    allocate (knots(nknots), coefs_bern(ncoefs))

    call ppolyfit (pp_bern, xx, ydat, deg, knots, coefs_bern, status)

    deallocate (xx)

    ncoefs = ppoly_get_ncoefs (pp, n=nknots, k=deg)
    allocate (coefs(ncoefs))

    ! Convert to power-basis
    call ppoly_transform_basis (pp_bern, coefs_bern, pp, coefs, status)

    ! Evaluate at original points
    allocate (yy(nknots))
    call ppolyval (pp, knots, coefs, knots, yy, status=status)

    is_ok = all_close (yy, ydat(1,:), atol=1.0e-10_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. is_ok, &
        'Eval at original points')

    deallocate (yy)

    ! Compare to points evaluated via bernstein basis at points other
    ! than original
    n = 101
    allocate (xx(n), yy(n), yy_bern(n))

    call linspace (xx, lb, ub)
    call ppolyval (pp_bern, knots, coefs_bern, xx, yy_bern, status=status)
    call ppolyval (pp, knots, coefs, xx, yy, status=status)

    is_ok = all_close (yy, yy_bern, atol=1.0e-10_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. is_ok, &
        'Eval at non-original points, comparing to Bernstein basis')

    deallocate (xx, yy, yy_bern)

    ! Compare 1st derivatives
    ncoefs_bern = ppoly_get_ncoefs (pp_bern_d1, n=nknots, k=deg-1)
    ncoefs = ppoly_get_ncoefs (pp_d1, n=nknots, k=deg-1)
    allocate (coefs_d1(ncoefs), coefs_bern_d1(ncoefs_bern))

    call ppolyder (pp_bern, knots, coefs_bern, pp_bern_d1, coefs_bern_d1, m=1, status=status)
    call ppolyder (pp, knots, coefs, pp_d1, coefs_d1, m=1, status=status)

    ! Compare values to bernstein values
    n = 201
    allocate (xx(n), yy(n), yy_bern(n))
    call linspace (xx, lb, ub)

    call ppolyval (pp_bern_d1, knots, coefs_bern_d1, xx, yy_bern, status=status)
    call ppolyval (pp_d1, knots, coefs_d1, xx, yy, status=status)

    is_ok = all_close (yy, yy_bern, atol=1.0e-10_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. is_ok, &
        'D1: Comparing to Bernstein value at non-original points')

    deallocate (xx, yy, yy_bern)

    ! Compute first deratives to true values at original knots as these
    ! were used to fit the data.
    allocate (yy(nknots))
    call ppolyval (pp_d1, knots, coefs_d1, knots, yy, status=status)
    is_ok = all_close (yy, ydat(2,:), atol=1.0e-8_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. is_ok, &
        'D1: Comparing to true values at original knots')

end subroutine

end
