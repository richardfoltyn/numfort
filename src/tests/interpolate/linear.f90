program test_interp_linear

    use iso_fortran_env
    use corelib_testing
    use numfort_arrays
    use numfort_interpolate
    implicit none

    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_interpolate unit tests")

    call test_exact (tests)
    call test_extrap (tests)
    call test_truncate (tests)
    call test_random (tests)

    ! print test statistics
    call tests%print ()

end subroutine

!>  Tests for linear extrapolation when x-coordinate of interpolated value
!   is outside of data points boundaries.
subroutine test_extrap (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(10) :: xp, fp, fx, x
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
        call interp_linear (xp(i), xp(j:k), fp(j:k), fx(i), ext=.true.)
    end do

    call tc%assert_true (all(abs(fx-fp) < 1d-12), &
        "Extrapolation, scalar argument")

    ! ===== Array arguments =====
    fx = 0.0
    call interp_linear (xp, xp(j:k), fp(j:k), fx, ext=.true.)
    call tc%assert_true (all(abs(fx-fp) < 1d-12), &
        "Extrapolation, array argument")

end subroutine

!>  Tests for "interpolating" exactly at function data points.
subroutine test_exact (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(10) :: xp, fp, fx, x
    integer :: i, j, k

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

    real (PREC), dimension(10) :: xp, fp, fx, x, fx1
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
        call interp_linear (xp(i), xp(j:k), fp(j:k), fx(i), ext=.false.)
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
        call interp_linear (xp(i), xp(j:k), fp(j:k), fx(i), ext=.false., left=left)
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
        call interp_linear (xp(i), xp(j:k), fp(j:k), fx(i), ext=.false., right=right)
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

    call interp_linear (xp, xp(j:k), fp(j:k), fx, ext=.false.)

    ! Compute expected result
    fx1 = fp
    call truncate (xp, xp(j), xp(k), fx1, fp(j), fp(k))

    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Truncation: array argument, default left/right values")

    ! >>> Truncation, exlicitly specify left value
    fx = 0.0
    left = -1.0

    call interp_linear (xp, xp(j:k), fp(j:k), fx, left=left)

    ! compute expected result
    fx1 = fp
    call truncate (xp, xp(j), xp(k), fx1, left, fp(k))

    call tc%assert_true (all(abs(fx-fx1) < 1d-12), &
        "Truncation: array argument, explicit left value")

    ! >>> Truncation: exlicitly specify right value
    fx = 0.0
    right = -10.0
    call interp_linear (xp, xp(j:k), fp(j:k), fx, right=right)

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

end program
