program test_interp_linear

    use, intrinsic :: iso_fortran_env
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

    call tests%set_label ("numfort_interpolate unit tests using OpenMP")

    call test_extrap (tests)

    ! print test statistics
    call tests%print ()

end subroutine

!>  Tests for linear extrapolation when x-coordinate of interpolated value
!   is outside of data points boundaries.
subroutine test_extrap (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: n = int(1e5)
    real (PREC), dimension(10) :: xp, fp
    real (PREC), dimension(:), allocatable :: fx_hat, x, fx
    integer :: i

    tc => tests%add_test ("Linear extrapolation")

    ! function f(x) = 2x + 1 on [10, 20]
    call linspace (xp, 1.0d1, 2.0d1)
    call func1 (xp, fp)

    allocate (x(n), fx(n), fx_hat(n))
    ! evaluate on [0,100], ie using extrapolation
    call linspace (x, 0.0d0, 1.0d2)
    fx = 0.0

    !$omp parallel default(shared) private(i)
    !$omp do
    do i = 1, size(x)
        call interp_linear (x(i), xp, fp, fx_hat(i), ext=.true.)
    end do
    !$omp end do

    !$omp end parallel

    ! compute true values f(x)
    call func1 (x, fx)

    call tc%assert_true (all(abs(fx-fx_hat) < 1d-12), &
        "Extrapolation, scalar argument")

end subroutine

pure subroutine func1 (x, fx)
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: fx

    fx = 2.0d0 * x + 1.0d0
end subroutine

end program
