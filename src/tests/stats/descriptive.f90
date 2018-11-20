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

    call tests%set_label ("numfort_stats_core unit tests")

    call test_degenerate (tests)
    call test_1d (tests)
    call test_2d (tests)
    call test_cov (tests)

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

    tc => tests%add_test ("Mean/SD degenerate inputs")

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



subroutine test_cov (tests)
    !*  Unit tests for COV and CORRCOEF routines
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(:,:), allocatable :: x, vcv, vcv_desired
    real (PREC), dimension(:), allocatable :: eps, var_x
    integer :: nvars, nobs
    real (PREC) :: s, alpha, beta1, beta2, var_eps1, var_eps2, diff
    type (status_t) :: status
    logical :: vcv_ok

    tc => tests%add_test ('Unit tests for COV and CORRCOEF routines')

    call set_seed (1234)

    ! === Input validation tests ===
    nobs = 1
    nvars = 1
    allocate (x(nobs, nvars), vcv(nvars, nvars+1))
    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid input: non-conformable X and COV')

    deallocate (x, vcv)

    ! Invalid DIM argument
    nobs = 1
    nvars = 1
    allocate (x(nobs,nvars), vcv(nobs,nobs))

    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, dim=0, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid input: DIM=0')

    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, dim=3, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid input: DIM=3')

    ! Invalid DOF argument
    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, dof=-1, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid input: DOF=0')

    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, dof=2, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid input: DOF=2')

    deallocate (x, vcv)

    ! === Degenerate input data ===
    nobs = 1
    nvars = 1
    allocate (x(nobs,nvars), vcv(nvars,nvars))
    x(:,:) = 1.0
    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. all(vcv == 0.0), &
        'COV called with a single observation')
    deallocate (x, vcv)

    ! Single input variable, multiple observations
    nobs = 10
    nvars = 1
    allocate (x(nobs,nvars), vcv(nvars,nvars))
    call random_number (x)
    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, status=status)
    call std (x(:,1), s)
    vcv_ok = all(abs(vcv - s**2.0_PREC) < 1.0e-8_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. vcv_ok, &
        'COV called with single variable, multiple observations')
    deallocate (x, vcv)

    ! Single input variable, swapped dimensions
    nobs = 2
    nvars = 1
    allocate (x(nvars,nobs), vcv(nvars,nvars))
    call random_number(x)
    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, dim=2, status=status)
    call std (x(1,:), s)
    vcv_ok = all(abs(vcv - s**2.0_PREC) < 1.0e-8_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. vcv_ok, &
        'COV called with single variable, DIM=2')
    deallocate (x, vcv)

    ! === Input data with many variables/obs ===
    ! v2 is exact linear function of v1
    nobs = 100
    nvars = 2
    allocate (x(nobs,nvars), vcv(nvars,nvars))
    call random_number (x(:,1))
    alpha = 0.123d0
    beta1 = 0.456d0
    x(:,2) = alpha + beta1 * x(:,1)
    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, status=status)
    ! Correct value of Cov(x1,x2) = Cov(x1,alpha+beta1*x1) = beta1*Var(x1)
    vcv_ok = abs(beta1*vcv(1,1) - vcv(1,2)) < 1.0e-8_PREC
    call tc%assert_true (status == NF_STATUS_OK .and. vcv_ok, &
        'COV called with X2 = alpha + beta1*X1')
    deallocate (x, vcv)

    ! Add some noise such that
    ! X2 = alpha + beta1 * X1 + eps
    nobs = 10000
    nvars = 2
    allocate (x(nvars,nobs), vcv(nvars,nvars), eps(nobs))
    call random_number (x(1,:))
    call rvs (norm, eps, scale=0.1_PREC)
    alpha = -0.2345d0
    beta1 = 1.234d0
    x(2,:) = alpha + beta1 * x(1,:) + eps
    status = NF_STATUS_UNDEFINED
    call cov (x, vcv, dim=2, status=status)
    ! Correct value of Cov(x1,x2) = Cov(x1,alpha+beta1*x1+eps) = beta1*Var(x1)
    ! as eps and x1 are uncorrelated
    vcv_ok = abs(beta1*vcv(1,1) - vcv(1,2)) < 1.0e-3_PREC
    call tc%assert_true (status == NF_STATUS_OK .and. vcv_ok, &
        'COV called with X2 = alpha + beta1*X1 + eps, DIM=2')
    deallocate (x, vcv, eps)

    ! === CORRCOEF tests ===
    ! Degenerate input data
    nobs = 1
    nvars = 1
    allocate (x(nobs, nvars), vcv(nvars,nvars))
    x(:,:) = 1.0
    status = NF_STATUS_UNDEFINED
    call corrcoef (x, vcv, status=status)
    vcv_ok = all(vcv == 1.0)
    call tc%assert_true (status == NF_STATUS_OK .and. vcv_ok, &
        'CORRCOEF called with 1 variable / 1 obs.')
    deallocate (x, vcv)

    ! Call with one variable, multiple obs.
    nobs = 21343
    nvars = 1
    allocate (x(nobs, nvars), vcv(nvars,nvars))
    call random_number (x)
    status = NF_STATUS_UNDEFINED
    call corrcoef (x, vcv, status=status)
    vcv_ok = all(vcv == 1.0)
    call tc%assert_true (status == NF_STATUS_OK .and. vcv_ok, &
        'CORRCOEF called with 1 variable, multiple obs.')
    deallocate (x, vcv)

    ! === CORRCOEF tests with several variables
    nobs = 50000
    nvars = 3
    allocate (x(nobs,nvars), vcv(nvars,nvars), eps(nobs))
    allocate (var_x(nvars))
    call random_number (x(:,1))
    call rvs (norm, eps, scale=0.1_PREC)
    call std (eps, var_eps1)
    var_eps1 = var_eps1**2.0

    beta1 = 0.234d0
    beta2 = -7.2389d0
    x(:,2) = 0.213d0 + beta1 * x(:,1) + eps

    call rvs (norm, eps, scale=0.05_PREC)
    call std (eps, var_eps2)
    var_eps2 = var_eps2**2.0

    x(:,3) = -12.234d0 + beta2 * x(:,1) + eps
    call std (x, var_x)
    var_x(:) = var_x**2.0

    status = NF_STATUS_UNDEFINED
    call corrcoef (x, vcv, status=status)
    ! Cov(x1,x2) = beta1 * var(x1) as above
    ! Var(x2) = beta1^2 * Var(x1) + Var(eps)
    ! Corr(x1,x2) = beta1 * Var(x1) / sqrt(Var(x1) * (beta^2*Var(x1) + Var(eps))
    ! For correlations between X2, X3:
    !   Cov(x2,x3)  = Cov(alpha1+beta1*x1+eps1, alpha2+beta2*x1+eps2) =
    !               = beta1 * beta2 * Var(x1)
    !   Corr(x2,x3) = (beta1*beta2*Var(X1))/sqrt((beta1^2*Var(x1)+Var(eps1)) *
    !                   (beta2^2*Var(X1)+Var(eps2))
    allocate (vcv_desired(nvars,nvars), source=1.0_PREC)
    vcv_desired(1,3) = beta2 * var_x(1) / sqrt(var_x(1) * (beta2**2.0*var_x(1) + var_eps2))
    vcv_desired(3,1) = vcv_desired(1,3)
    vcv_desired(2,1) = beta1 * var_x(1) / sqrt(var_x(1) * (beta1**2.0*var_x(1) + var_eps1))
    vcv_desired(1,2) = vcv_desired(2,1)
    vcv_desired(2,3) = beta1*beta2*var_x(1)/sqrt((beta1**2.0*var_x(1) + var_eps1)) &
        / sqrt(beta2**2.0 * var_x(1) + var_eps2)
    vcv_desired(3,2) = vcv_desired(2,3)

    diff = maxval(abs(vcv - vcv_desired))

    vcv_ok = diff < 1.0e-3_PREC
    call tc%assert_true (status == NF_STATUS_OK .and. vcv_ok, &
        'CORRCOEF called with X2, X3 linear functions of X1 with error terms')


end subroutine



end program
