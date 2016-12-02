program test_numfort_stats_dnorm

    use iso_fortran_env

    use corelib_strings
    use corelib_testing
    use numfort_arrays
    use numfort_stats, only: dnorm, norm

    implicit none

    real (real64), parameter :: tol = 1d-12

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_stats_dnorm unit tests")

    call test_cdf (tests)
    ! call test_sub2ind (tests)

    ! print test statistics
    call tests%print ()

end subroutine


subroutine test_cdf (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    type (str) :: msg
    real (real64), dimension(:), allocatable :: x, fx
    real (real64) :: x1, fx1
    real (real64), parameter, dimension(4) :: &
        means = [0.0d0, -21.0d0, 1d-5, 1.234d4], &
        sd = [1.0d0, 1d-4, 1d1, 1.2345d4]

    integer :: i, n

    tc => tests%add_test ("dnorm CDF test cases")

    do i = 1, size(means)

        msg = " (mean=" // str(means(i), 'g8.1e2') // ", sd=" // &
            str(sd(i), 'g8.3e2') // ")"

        x1 = means(i)
        call norm%cdf (x1, fx1, mean=means(i), sd=sd(i))
        call tc%assert_true (abs(fx1 - 0.5d0) < tol, &
            "Check CDF(mean) == 0.5" // msg)

        ! Check boundary behavior
        x1 = huge(x1)
        call norm%cdf (x1, fx1, mean=means(i), sd=sd(i))
        call tc%assert_true (abs(fx1 - 1.0d0) < tol, &
            "Check CDF(inf) == 1.0" // msg)

        x1 = -huge(x1)
        call norm%cdf (x1, fx1, mean=means(i), sd=sd(i))
        call tc%assert_true (abs(fx1) < tol, &
            "Check CDF(-inf) = 0.0" // msg)

        ! Test that CDF is non-decreasing
        n = 100
        allocate (x(n), fx(n))

        call linspace (x, -1d10, 1d10)
        call norm%cdf (x, fx, mean=means(i), sd=sd(i))
        call tc%assert_true (all(fx(2:)-fx(1:n-1) >= 0.0d0), &
            "Check that CDF is non-decreasing" // msg)

        deallocate (x, fx)
    end do

end subroutine


end program
