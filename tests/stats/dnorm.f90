program test_numfort_stats_dnorm

    use iso_fortran_env

    ! avoid nameclashes with CORELIB status_t type
    use fcore_common, FC_status_t => status_t
    use fcore_testing
    use numfort_arrays
    use numfort_common
    use numfort_stats, only: dnorm => dnorm_real64, norm, mean, std, cdf, pdf, rvs

    implicit none

    integer, parameter :: PREC = real64

    real (PREC), parameter :: tol = 1d-12

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_stats_dnorm unit tests")

    call test_cdf (tests)
    call test_cdf_rvs (tests)

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

        msg = " (loc=" // str(means(i), 'g8.1e2') // ", scale=" // &
            str(sd(i), 'g8.3e2') // ")"

        x1 = means(i)
        fx1 = cdf (norm, x1, loc=means(i), scale=sd(i))
        call tc%assert_true (abs(fx1 - 0.5d0) < tol, &
            "Check CDF(mean) == 0.5" // msg)

        ! Check boundary behavior
        x1 = huge(x1)
        fx1 = cdf (norm, x1, loc=means(i), scale=sd(i))
        call tc%assert_true (abs(fx1 - 1.0d0) < tol, &
            "Check CDF(inf) == 1.0" // msg)

        x1 = -huge(x1)
        fx1 = cdf (norm, x1, loc=means(i), scale=sd(i))
        call tc%assert_true (abs(fx1) < tol, &
            "Check CDF(-inf) = 0.0" // msg)

        ! Test that CDF is non-decreasing
        n = 100
        allocate (x(n), fx(n))

        call linspace (x, -1d10, 1d10)
        fx = cdf (norm, x, loc=means(i), scale=sd(i))
        call tc%assert_true (all(fx(2:)-fx(1:n-1) >= 0.0d0), &
            "Check that CDF is non-decreasing" // msg)

        deallocate (x, fx)
    end do

end subroutine

subroutine test_cdf_rvs (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    type (str) :: msg
    real (PREC), dimension(:), allocatable :: x, fx
    real (PREC), parameter, dimension(4) :: &
        means = [real (PREC) :: 0.0d0, -21.0d0, 1d-5, 1.234d4], &
        sd = [real (PREC) :: 1.0d0, 1d-4, 1d1, 1.2345d4]
    real (PREC) :: m, s, s2
    integer :: i, n
    type (status_t) :: status
    logical :: in_range

    tc => tests%add_test ("dnorm CDF(random) test cases")

    n = 100000
    allocate (x(n), fx(n))
    do i = 1, size(means)

        msg = " (loc=" // str(means(i), 'g8.1e2') // ", scale=" // &
            str(sd(i), 'g8.3e2') // ")"

        ! draw random sample of normally distributed RV
        call rvs (norm, x, loc=means(i), scale=sd(i))
        ! compute CDF(x): is standard uniformly distributed RV
        fx = cdf (norm, x, loc=means(i), scale=sd(i))

        in_range = all(fx >= 0.0_PREC) .and. all(fx <= 1.0_PREC)

        ! compute mean and variance
        call std (fx, s, m, status=status)
        s2 = s ** 2
        call tc%assert_true (in_range .and. abs(m-0.5) < 1d-3 .and. abs(s2 - 1.0d0/12) < 1d-3, &
            "Check that CDF(x) close to std. uniform" // msg)

    end do

    deallocate (x, fx)

end subroutine

end program
