program test_nf_stats_desc

    use iso_fortran_env

    use corelib_strings, only: str
    use corelib_testing, only: test_suite, test_case
    use numfort_arrays
    use numfort_stats
    use numfort_common

    implicit none

    integer, parameter :: PREC = real64
    real (real64), parameter :: tol = 1d-12

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_stats_dnorm unit tests")

    call test_degenerate (tests)
    call test_1d (tests)
    call test_2d (tests)

    ! print test statistics
    call tests%print ()

end subroutine


subroutine test_degenerate (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC) :: x(0), x10(1,0), x01(0,1), x1(1), x21(2,1), x12(1,2)
    real (PREC) :: m, m1(2)
    real (PREC) :: s, s1(2)
    type (status_t) :: status

    tc => tests%add_test ("Mean/SD degenerator inputs")

    call mean (x, m, status=status)
    call tc%assert_true (NF_STATUS_INVALID_ARG .in. status, &
        "Degenerate 1d input array")

    call mean (x10, m1, status=status)
    call tc%assert_true ((NF_STATUS_INVALID_ARG .in. status), &
        "Degenerate 2d input array (dim=2)")

    call mean (x01, m1, status=status)
    call tc%assert_true ((NF_STATUS_INVALID_ARG .in. status), &
        "Degenerate 2d input array (dim=1)")

    x1 = 1.234
    call mean (x1, m, status=status)
    call tc%assert_true (m == x1(1) .and. (NF_STATUS_OK .in. status), &
        "mean: 1d input, size(x) = 1")

    m = -1
    call std (x1, s, m, status=status)
    call tc%assert_true (m == x1(1) .and. s == 0.0_PREC .and. (NF_STATUS_OK .in. status), &
        "std: 1d input, size(x) = 1")


    x21(:, 1) = [7.234, -0.223423]
    call mean (x21, m1, dim=2, status=status)
    call tc%assert_true (all(x21(:, 1) == m1) .and. (NF_STATUS_OK .in. status), &
        "mean: 2d input of shape [2, 1]")

    x12(1, :) = [789.213, -2.123]
    call mean (x12, m1, dim=1, status=status)
    call tc%assert_true (all(x12(1, :) == m1) .and. (NF_STATUS_OK .in. status), &
        "mean: 2d input of shape [1, 2]")

    x21(:, 1) = [1.2345, -0.3214]
    s1 = -1.0
    call std(x21, s1, m1, dim=2, status=status)
    call tc%assert_true (all(x21(:, 1) == m1) .and. all(s1 == 0.0_PREC) .and. (NF_STATUS_OK .in. status), &
        "std: 2d input of shape [2, 1]")

    s1 = -1.0
    x12(1, :) = [0.213, -234.123]
    call std(x12, s1, m1, dim=1, status=status)
    call tc%assert_true (all(x12(1, :) == m1) .and. all(s1 == 0.0_PREC) .and. (NF_STATUS_OK .in. status), &
        "std: 2d input of shape [1, 2]")

end subroutine

subroutine test_1d (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: N = 101
    real (PREC) :: x(N)
    real (PREC) :: m, m2, s, s2
    integer :: i
    type (status_t) :: status

    tc => tests%add_test ("1d input arrays")

    x = [real (PREC) :: (i, i = 1, N)]
    call mean (x, m, status=status)
    call tc%assert_true (abs(m - (N+1)/2) < 1d-15 .and. (NF_STATUS_OK .in. status), &
        "Mean of 1d sequence")

    ! compute std. deviation manually (numerical instability should be
    ! negligible in this case)
    s2 = sqrt(sum((x -m) ** 2) / (N-1))
    call std (x, s, m2, status=status)
    call tc%assert_true (m == m2 .and. abs(s-s2) < 1d-15 .and. (NF_STATUS_OK .in. status), &
        "Std of 1d sequence")

end subroutine

subroutine test_2d (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: N = 100
    real (PREC), dimension(N,N) :: x
    real (PREC), dimension(N) :: m, m2, s, s2
    integer :: i
    type (status_t) :: status

    tc => tests%add_test ("2d input arrays")

    ! draw some random sample
    call random_number (x)

    ! test reduction over dim=1
    do i = 1, N
        m2(i) = sum(x(:, i)) / N
        s2(i) = sqrt(sum((x(:, i) - m2(i)) ** 2) / (N-1))
    end do

    ! Compute using library routines
    call mean (x, m, dim=1, status=status)
    call tc%assert_true (all(abs(m-m2) < tol) .and. (NF_STATUS_OK .in. status), &
        "Mean of 2d array along dim=1")

    call std (x, s, m, dim=1, status=status)
    call tc%assert_true (all(abs(m-m2) < tol) .and. all(abs(s-s2) < tol) .and. (NF_STATUS_OK .in. status), &
        "Std of 2d array along dim=1")

    ! test reduction over dim=2
    do i = 1, N
        m2(i) = sum(x(i, :)) / N
        s2(i) = sqrt(sum((x(i, :) - m2(i)) ** 2) / (N-1))
    end do

    ! Compute using library routines
    call mean (x, m, dim=2, status=status)
    call tc%assert_true (all(abs(m-m2) < tol) .and. (NF_STATUS_OK .in. status), &
        "Mean of 2d array along dim=2")

    call std (x, s, m, dim=2, status=status)
    call tc%assert_true (all(abs(m-m2) < tol) .and. all(abs(s-s2) < tol) .and. (NF_STATUS_OK .in. status), &
        "Std of 2d array along dim=2")
end subroutine

end program
