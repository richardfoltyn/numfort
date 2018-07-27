program test_nf_stats_desc

    use iso_fortran_env

    use fcore_testing, only: test_suite, test_case
    use numfort_arrays
    use numfort_stats
    use numfort_common
    use numfort_common_testing
    use numfort_interpolate

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

    call test_quantile (tests)

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



subroutine test_quantile (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC) :: wgt_lb
    real (PREC), dimension(:), allocatable :: bins, pmf, cdf
    real (PREC), dimension(:), allocatable :: rnk, q, q_ok
    integer :: n, nq, i
    type (status_t) :: status

    tc => tests%add_test ("Unit tests for PERCENTILE routine")

    call set_seed (1234)

    ! ==== Test 1 =====
    n = 5
    allocate (bins(n+1), pmf(n))
    call linspace (bins, 0.0_PREC, real(n, PREC))
    pmf(:) = 1.0_PREC / n

    allocate (cdf(n+1))
    cdf(1) = 0.0
    do i = 1, n
        cdf(i+1) = cdf(i) + pmf(i)
    end do

    nq = 5
    allocate (rnk(nq), q(nq), q_ok(nq))
    ! Test at "interior" quantiles that do not coindice with inverse CDF
    ! values on bin edges
    wgt_lb = 0.25d0
    rnk(:) = (1.0-wgt_lb) * cdf(2:n+1) + wgt_lb * cdf(1:n)

    ! Linear interpolation
    q_ok(:) = wgt_lb * bins(1:n) + (1.0-wgt_lb) * bins(2:n+1)

    call quantile (bins, pmf, rnk, q, interp='linear', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 1: interp='linear'")

    ! Lower "interpolation"
    q_ok = bins(1:n)
    call quantile (bins, pmf, rnk, q, interp='lower', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 1: interp='lower'")

    ! Higher "interpolation"
    q_ok = bins(2:n+1)
    call quantile (bins, pmf, rnk, q, interp='higher', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 1: interp='higher'")

    ! Midpoint interpolation
    wgt_lb = 0.5d0
    q_ok = wgt_lb * bins(1:n) + (1.0-wgt_lb) * bins(2:n+1)
    call quantile (bins, pmf, rnk, q, interp='midpoint', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 1: interp='midpoint'")

    ! ===== Test 2 =====
    ! Test with 0 elements in PMF at the extremes
    pmf(:) = 0.0
    pmf(2:4) = 1.0_PREC / 3.0

    ! Lower interpolation where 0, 1.0 are in RNK
    rnk(:) = [0.0d0, 1.0d0/6.0d0, 3.0d0/6.0d0, 5.0d0/6.0d0, 1.0d0]
    q_ok(1:2) = bins(2)
    q_ok(3) =  bins(3)
    q_ok(4) = bins(4)
    q_ok(5) = bins(5)

    call quantile (bins, pmf, rnk, q, interp='lower', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 2: interp='lower', PMF with zeros")

    ! Higher interpolation where 0.0, 1.0 are in RNK
    q_ok(1) = bins(2)
    q_ok(2) = bins(3)
    q_ok(3) = bins(4)
    q_ok(4) = bins(5)
    q_ok(5) = bins(5)

    call quantile (bins, pmf, rnk, q, interp='higher', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 2: interp='higher', PMF with zeros")

    ! Linear interpolation
    q_ok(1) = bins(2)
    q_ok(2:4) = (bins(2:4) + bins(3:5)) / 2.0
    q_ok(5) = bins(5)

    call quantile (bins, pmf, rnk, q, interp='linear', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 2: interp='linear', PMF with zeros")

    ! Midpoint
    q_ok(1) = q_ok(2)
    q_ok(5) = q_ok(4)
    call quantile (bins, pmf, rnk, q, interp='midpoint', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 2: interp='midpoint', PMF with zeros")

    ! ===== Test 3 ======
    ! PMF with zeros in middle
    pmf(:) = 1.0
    pmf(3:4) = 0.0
    pmf(:) = pmf / sum(pmf)

    rnk(:) = [0.0d0, 1.0d0/6.0d0, 3.0d0/6.0d0, 5.0d0/6.0d0, 1.0d0]

    ! Linear interpolation
    q_ok(1) = bins(1)
    q_ok(2:3) = (bins(1:2) + bins(2:3)) / 2.0
    q_ok(4) = (bins(5) + bins(6)) / 2.0
    q_ok(5) = bins(6)
    call quantile (bins, pmf, rnk, q, interp='linear', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 3: interp='linear', PMF with interior zeros")

    ! lower "interpolation"
    q_ok(1:2) = bins(1)
    q_ok(3) = bins(2)
    q_ok(4) = bins(5)
    q_ok(5) = bins(6)

    call quantile (bins, pmf, rnk, q, interp='lower', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 3: interp='lower', PMF with interior zeros")

    ! higher "interpolation"
    q_ok(1) = bins(1)
    q_ok(2) = bins(2)
    q_ok(3) = bins(3)
    q_ok(4:5) = bins(6)
    call quantile (bins, pmf, rnk, q, interp='higher', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 3: interp='higher', PMF with interior zeros")

    ! midpoint interpolation
    q_ok(1:2) = (bins(1) + bins(2)) / 2.0
    q_ok(3) = (bins(2) + bins(3)) / 2.0
    q_ok(4:5) = (bins(5) + bins(6)) / 2.0
    call quantile (bins, pmf, rnk, q, interp='midpoint', status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (q, q_ok, atol=1.0d-10), &
        "Test 3: interp='midpoint', PMF with interior zeros")


end subroutine

end program
