program test_interp_linear
    !*  Unit tests for linear interpolation in various dimensions.

    use iso_fortran_env
    use numfort_arrays
    use numfort_common_testing
    use numfort_interpolate

    use fcore_testing

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("LINEAR interpolation unit tests")

    call test_exact (tests)
    call test_extrap (tests)
    call test_truncate (tests)
    call test_random (tests)

    call test_bilinear_input_checks (tests)
    call test_bilinear_interpolate (tests)
    call test_bilinear_extrapolate (tests)

    ! print test statistics
    call tests%print ()

end subroutine

!>  Tests for linear extrapolation when x-coordinate of interpolated value
!   is outside of data points boundaries.
subroutine test_extrap (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(10) :: xp, fp, fx
    integer :: i, j, k

    tc => tests%add_test ("Linear extrapolation")

    ! function f(x) = 2x + 1 on [0, 1]
    call linspace (xp, 0.0d0, 1.0d0)
    fp = 2 * xp + 1.0d0

    ! ===== Scalar arguments =====
    ! evaluate outside of domain with extrapolation
    fx = 0.0
    ! limit xp/fp to subrange of values
    j = size(xp) / 2 - 1
    k = size(xp) / 2 + 1

    do i = 1, size(xp)
        call interp_linear (xp(i), xp(j:k), fp(j:k), fx(i), &
            ext=NF_INTERP_EVAL_EXTRAPOLATE)
    end do

    call tc%assert_true (all(abs(fx-fp) < 1d-12), &
        "Extrapolation, scalar argument")

    ! ===== Array arguments =====
    fx = 0.0
    call interp_linear (xp, xp(j:k), fp(j:k), fx, ext=NF_INTERP_EVAL_EXTRAPOLATE)
    call tc%assert_true (all(abs(fx-fp) < 1d-12), &
        "Extrapolation, array argument")

end subroutine

!>  Tests for "interpolating" exactly at function data points.
subroutine test_exact (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(10) :: xp, fp, fx, x
    integer :: i

    tc => tests%add_test ("Exact 'interpolation': x in xp")

    ! function f(x) = 2x + 1 on [0, 1]
    call linspace (xp, 0.0d0, 1.0d0)
    fp = 2 * xp + 1.0d0

    ! ===== Scalar arguments =====
    ! evaluate at all points in xp
    fx = 0.0
    do i = 1, size(xp)
        call interp_linear (xp(i), xp, fp, fx(i))
    end do
    call tc%assert_true (all(abs(fx-fp) < 1d-12), &
        "Exact interpolation, scalar argument")

    ! ===== Array arguments =====
    x = xp
    fx = 0.0d0

    call interp_linear (x, xp, fp, fx)
    call tc%assert_true (all(abs(fp-fx) < 1d-12), &
        "Exact interpolation, array argument")

end subroutine

!>  Tests for interpolating at random points within the boundaries defined
!   but data points.
subroutine test_random (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(10) :: xp, fp, fx, x, fx1
    integer :: i
    real (PREC) :: x1, xn

    x1 = 0.0d0
    xn = 1.0d0

    tc => tests%add_test ("Linear function f(x) with random x")

    ! function f(x) = 2x + 1 on [0, 1]
    call linspace (xp, x1, xn)
    fp = 2 * xp + 1.0d0

    call random_number (x)
    ! correct function values
    x = x * (xn-x1) + x1
    fx1 = 2 * x + 1.0d0

    ! ===== Scalar arguments =====
    ! evaluate at all points in xp
    fx = 0.0
    do i = 1, size(xp)
        call interp_linear (x(i), xp, fp, fx(i))
    end do

    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Scalar argument")

    ! ===== Array arguments =====
    fx = 0.0d0

    call interp_linear (x, xp, fp, fx)
    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Array argument")

end subroutine

!>  Tests for interpolation at x-coordinates outside of data point boundaries
!   when linear extrapolation is disabled.
subroutine test_truncate (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(10) :: xp, fp, fx, fx1
    integer :: i, j, k
    real (PREC) :: left, right

    tc => tests%add_test ("Truncation: left/right values for x outside of boundaries")

    ! function f(x) = 2x + 1 on [0, 1]
    call linspace (xp, 0.0d0, 1.0d0)
    fp = 2 * xp + 1.0d0

    ! ===== Scalar arguments =====

    ! >>> Truncate outside of specified domain, default left/right values
    fx = 0.0
    ! limit xp/fp to subrange of values
    j = size(xp) / 2 - 1
    k = size(xp) / 2 + 1

    do i = 1, size(xp)
        call interp_linear (xp(i), xp(j:k), fp(j:k), fx(i), &
            ext=NF_INTERP_EVAL_BOUNDARY)
    end do

    ! Compute expected result
    fx1 = fp
    call truncate (xp, xp(j), xp(k), fx1, fp(j), fp(k))

    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Truncation: scalar argument, default left/right values")

    ! >>> Truncation, exlicitly specify left value
    fx = 0.0
    left = -1.0
    do i = 1, size(xp)
        call interp_linear (xp(i), xp(j:k), fp(j:k), fx(i), &
            ext=NF_INTERP_EVAL_CONST, left=left)
    end do

    ! compute expected result
    fx1 = fp
    call truncate (xp, xp(j), xp(k), fx1, left, fp(k))

    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Truncation: scalar argument, explicit left value")

    ! >>> Truncation: exlicitly specify right value
    fx = 0.0
    right = -10.0
    do i = 1, size(xp)
        call interp_linear (xp(i), xp(j:k), fp(j:k), fx(i), &
            ext=NF_INTERP_EVAL_CONST, right=right)
    end do

    ! compute expected result
    fx1 = fp
    call truncate (xp, xp(j), xp(k), fx1, fp(j), right)

    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Truncation: scalar argument, explicit right value")


    ! ===== Array arguments =====

    ! >>> Truncate outside of specified domain, default left/right values
    fx = 0.0
    ! limit xp/fp to subrange of values
    j = size(xp) / 2 - 1
    k = size(xp) / 2 + 1

    call interp_linear (xp, xp(j:k), fp(j:k), fx, ext=NF_INTERP_EVAL_BOUNDARY)

    ! Compute expected result
    fx1 = fp
    call truncate (xp, xp(j), xp(k), fx1, fp(j), fp(k))

    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Truncation: array argument, default left/right values")

    ! >>> Truncation, exlicitly specify left value
    fx = 0.0
    left = -1.0

    call interp_linear (xp, xp(j:k), fp(j:k), fx, ext=NF_INTERP_EVAL_CONST, &
        left=left)

    ! compute expected result
    fx1 = fp
    call truncate (xp, xp(j), xp(k), fx1, left, fp(k))

    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Truncation: array argument, explicit left value")

    ! >>> Truncation: exlicitly specify right value
    fx = 0.0
    right = -10.0
    call interp_linear (xp, xp(j:k), fp(j:k), fx, ext=NF_INTERP_EVAL_CONST, &
        right=right)

    ! compute expected result
    fx1 = fp
    call truncate (xp, xp(j), xp(k), fx1, fp(j), right)

    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Truncation: array argument, explicit right value")

end subroutine

subroutine truncate (x, lb, ub, fx, left, right)
    real (PREC), dimension(:) :: x, fx
    real (PREC), intent(in) :: lb, ub, left, right

    where (x < lb)
        fx = left
    else where (x  > ub)
        fx = right
    end where
end subroutine


subroutine test_bilinear_input_checks (tests)
    !*  Checks input verification for bilinear interpolation routines.
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: xp1, xp2, x1, x2, fx1d
    real (PREC) :: fx
    real (PREC), dimension(:,:), allocatable :: fp
    type (status_t) :: status

    tc => tests%add_test ("Bilinear interpolation: input checks")

    ! Test with grid arrays with less than 2 elements
    allocate (xp1(1), xp2(1), fp(1,1))
    status = NF_STATUS_OK
    call interp_bilinear (1.0_PREC, 1.0_PREC, xp1, xp2, fp, fx, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Grid point arrays have fewer than two elements")
    deallocate (xp1, xp2, fp)

    ! Test with incompatible grid arrays
    allocate (xp1(2), xp2(3), fp(2,2))
    status = NF_STATUS_OK
    call interp_bilinear (1.0_PREC, 1.0_PREC, xp1, xp2, fp, fx, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Grid point arrays are incompatible")
    deallocate (xp1, xp2, fp)

    ! Test with eval input arrays of different size
    allocate (xp1(2), xp2(2), fp(2,2))
    allocate (x1(1), x2(2), fx1d(1))
    status = NF_STATUS_OK
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx1d, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Eval input arrays are incompatible")
    deallocate (xp1, xp2, fp, x1, x2, fx1d)

    ! Test with eval output array of different size
    allocate (xp1(2), xp2(2), fp(2,2))
    allocate (x1(2), x2(2), fx1d(1))
    status = NF_STATUS_OK
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx1d, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Eval input arrays are incompatible")
    deallocate (xp1, xp2, fp, x1, x2, fx1d)

end subroutine



subroutine test_bilinear_interpolate (tests)
    !*  Tests for bilinear interpolation at interior points
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(:), allocatable :: xp1, xp2, x1, x2
    real (PREC), dimension(:,:), allocatable :: fp
    real (PREC), dimension(:), allocatable :: fx, fx_ok
    integer :: n1, n2, n, i, j
    logical :: res_ok
    type (status_t) :: status

    tc => tests%add_test ("Bilinear interpolation tests (interior)")

    n1 = 10
    n2 = 11
    allocate (xp1(n1), xp2(n2), fp(n1,n2))

    call linspace (xp1, 0.0_PREC, 10.0_PREC)
    call powerspace (xp2, 0.0_PREC, 5.0_PREC, 2)

    forall (j=1:size(xp2)) fp(:,j) = func1_2d (xp1, xp2(j))

    ! Test with 1d array interface
    n = 15
    allocate (x1(n), x2(n), fx_ok(n), fx(n))
    call random_number (x1)
    call random_number (x2)

    ! Rescale to approx. cover most parts of rectangular grid
    x1(:) = x1 * 10.0
    x2(:) = x2 * 5.0

    ! Correct function values
    fx_ok(:) = func1_2d(x1, x2)

    ! Interpolated values
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx, status=status)
    res_ok = (status == NF_STATUS_OK) .and. all_close (fx, fx_ok, atol=1.0e-10_PREC)
    call tc%assert_true (res_ok, "Interpolate at random interior points; array input")

    ! Test scalar interface
    fx(:) = 0.0
    fx_ok(1) = func1_2d (x1(1), x2(1))
    call interp_bilinear (x1(1), x2(1), xp1, xp2, fp, fx(1), status=status)
    res_ok = (status == NF_STATUS_OK) .and. all_close (fx(1), fx_ok(1), atol=1.0e-10_PREC)
    call tc%assert_true (res_ok, "Interpolate at random interior point; scalar input")

    deallocate (x1, x2, fx, fx_ok)

    ! Test evaluating exactly at grid points
    n = n1 * n2
    allocate (x1(n), x2(n), fx(n), fx_ok(n))
    do j = 1, size(xp2)
        i = (j-1) * size(xp1)
        x1(i+1:i+n1) = xp1
        x2(i+1:i+n1) = xp2(j)
    end do

    fx_ok(:) = func1_2d (x1, x2)

    call interp_bilinear (x1, x2, xp1, xp2, fp, fx, status=status)

    res_ok = .true.
    do j = 1, size(xp2)
        i = (j-1) * n1
        res_ok = res_ok .and. all_close (fx(i+1:i+n1), fp(:,j), atol=1.0e-10_PREC)
    end do
    call tc%assert_true (res_ok, "Interpolate exactly at grid points")


    deallocate (x1, x2, fx, fx_ok)
end subroutine



subroutine test_bilinear_extrapolate (tests)
    !*  Tests for bilinear interpolation when not all points are interior
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(:), allocatable :: xp1, xp2, x1, x2
    real (PREC), dimension(:,:), allocatable :: fp
    real (PREC), dimension(:), allocatable :: fx, fx_ok
    integer :: n1, n2, n, i, j
    logical :: res_ok
    type (status_t) :: status
    real (PREC) :: fill_value

    tc => tests%add_test ("Bilinear interpolation tests (non-interior)")

    n1 = 10
    n2 = 11
    allocate (xp1(n1), xp2(n2), fp(n1,n2))

    call linspace (xp1, -1.0_PREC, 1.0_PREC)
    call powerspace (xp2, -1.0_PREC, 1.0_PREC, 2)

    forall (j=1:size(xp2)) fp(:,j) = func1_2d (xp1, xp2(j))

    ! Test with 1d array interface
    n = 10
    allocate (x1(n), x2(n), fx_ok(n), fx(n))
    call random_number (x1)
    call random_number (x2)

    ! Rescale to approx. cover most parts of rectangular grid
    x1(:) = x1 * 2.0 - 1.0
    x2(:) = x2 * 2.0 - 1.0

    ! Set some non-interior points (both or either dimension)
    x1(1) = -1.5
    x2(2) = -1.1
    x1(3) = 1.234d0
    x2(4) = 1.456d0
    x1(5) = -1.234d0
    x2(5) = -2.345d0
    x1(6) = 2.3456d0
    x2(6) = 3.4456d0

    fill_value = -1000.0

    ! Correct function values
    fx_ok(:) = func1_2d (x1, x2)

    ! Overwrite non-interior points with FILL_VALUE
    do i = 1,n
        if (x1(i) < xp1(1) .or. x1(i) > xp1(n) .or. x2(i) < xp2(1) .or. x2(i) > xp2(n)) then
            fx_ok(i) = fill_value
        end if
    end do

    ! Interpolated values; set non-internal values to -1000
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx, fill_value=fill_value, status=status)
    res_ok = (status == NF_STATUS_OK) .and. all_close (fx, fx_ok, atol=1.0e-10_PREC)
    call tc%assert_true (res_ok, "Interpolation at non-int. values using FILL_VALUE")

    ! Test with extrapolation
    fx_ok(:) = func1_2d (x1, x2)
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx, status=status)
    res_ok = (status == NF_STATUS_OK) .and. all_close (fx, fx_ok, atol=1.0e-10_PREC)
    call tc%assert_true (res_ok, "Interpolation at non-int. values using extrapolation")

    deallocate (x1, x2, fx, fx_ok)

end subroutine




elemental function func1_2d(x1, x2) result(fx)
    !*  Function f: R^2 -> R that can be interpolated exactly by bilinear
    !   interpolation
    real (PREC), intent(in) :: x1
    real (PREC), intent(in) :: x2
    real (PREC) :: fx

    fx = 2.0 * x1 - 3.0 * x2 + 0.5d0 * x1 * x2 + 1.0
end function

end program
