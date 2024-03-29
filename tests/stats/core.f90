

program test_numfort_stats_core

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

    call test_mean_std_input_checks (tests)
    call test_mean_std_degenerate (tests)
    call test_mean_std_1d (tests)
    call test_mean_std_2d (tests)

    call test_standardize_1d (tests)
    call test_standardize_2d (tests)

    call test_cov (tests)

    call test_quantile (tests)

    ! print test statistics
    call tests%print ()

end subroutine


subroutine test_mean_std_degenerate (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC) :: x(0), x10(1,0), x01(0,1), x1(1), x21(2,1), x12(1,2)
    real (PREC) :: m, m1(2)
    real (PREC) :: s, s1(2)
    type (status_t) :: status

    tc => tests%add_test ("Mean/SD degenerate inputs")

    ! --- zero-obs. inputs ---

    status = NF_STATUS_UNDEFINED
    call mean (x, m, status=status)
    call tc%assert_true (NF_STATUS_INVALID_ARG .in. status, &
        "Degenerate 1d input array")

    status = NF_STATUS_UNDEFINED
    call mean (x10, m1(1:1), dim=2, status=status)
    call tc%assert_true ((NF_STATUS_INVALID_ARG .in. status), &
        "2d input array, Nvars > 0, Nobs = 0 (dim=2)")

    status = NF_STATUS_UNDEFINED
    call mean (x01, m1(1:1), dim=1, status=status)
    call tc%assert_true ((NF_STATUS_INVALID_ARG .in. status), &
        "2d input array, Nvars > 0, Nobs=0 (dim=1)")

    status = NF_STATUS_UNDEFINED
    call std (x, s, status=status)
    call tc%assert_true (NF_STATUS_INVALID_ARG .in. status, &
        "STD: degenerate 1d input array")

    status = NF_STATUS_UNDEFINED
    call std (x10, s1(1:1), dim=2, status=status)
    call tc%assert_true ((NF_STATUS_INVALID_ARG .in. status), &
        "2d input array, Nvars > 0, Nobs = 0 (dim=2)")

    status = NF_STATUS_UNDEFINED
    call std (x01, s1(1:1), dim=1, status=status)
    call tc%assert_true ((NF_STATUS_INVALID_ARG .in. status), &
        "2d input array, Nvars > 0, Nobs=0 (dim=1)")

    ! Note: passing zero variables but non-zero observations should return
    ! NF_STATUS_OK, since the problem is well-defined for all variables
    ! in this empty set!
    status = NF_STATUS_UNDEFINED
    call std (x10, s1(1:0), dim=1, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "2d input array, Nvars = 0, Nobs > 0 (dim=1)")

    ! Note: passing zero variables but non-zero observations should return
    ! NF_STATUS_OK, since the problem is well-defined for all variables
    ! in this empty set!
    status = NF_STATUS_UNDEFINED
    call std (x01, s1(1:0), dim=2, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "2d input array, Nvars = 0, Nobs > 0 (dim=2)")


    ! --- zero-variable inputs ---

    status = NF_STATUS_UNDEFINED
    call mean (x10, m1(1:0), dim=1, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "2d input array, Nvars = 0, Nobs > 0 (dim=1)")

    status = NF_STATUS_UNDEFINED
    call mean (x01, m1(1:0), dim=2, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "2d input array, Nvars = 0, Nobs > 0 (dim=2)")

    ! --- 1d input of length 1 ---

    x1 = 1.234
    status = NF_STATUS_UNDEFINED
    call mean (x1, m, status=status)
    call tc%assert_true (m == x1(1) .and. (NF_STATUS_OK .in. status), &
        "mean: 1d input, size(x) = 1")

    m = -1
    status = NF_STATUS_UNDEFINED
    call std (x1, s, m, status=status)
    call tc%assert_true (m == x1(1) .and. s == 0.0_PREC .and. (NF_STATUS_OK .in. status), &
        "std: 1d input, size(x) = 1")

    ! --- 2d input with one dim. of lenth 1 ---

    x21(:, 1) = [7.234, -0.223423]
    status = NF_STATUS_UNDEFINED
    call mean (x21, m1, dim=2, status=status)
    call tc%assert_true (all(x21(:, 1) == m1) .and. (NF_STATUS_OK .in. status), &
        "mean: 2d input of shape [2, 1]")

    x12(1, :) = [789.213, -2.123]
    status = NF_STATUS_UNDEFINED
    call mean (x12, m1, dim=1, status=status)
    call tc%assert_true (all(x12(1, :) == m1) .and. (NF_STATUS_OK .in. status), &
        "mean: 2d input of shape [1, 2]")

    x21(:, 1) = [1.2345, -0.3214]
    s1 = -1.0
    status = NF_STATUS_UNDEFINED
    call std(x21, s1, m1, dim=2, status=status)
    call tc%assert_true (all(x21(:, 1) == m1) .and. all(s1 == 0.0_PREC) .and. (NF_STATUS_OK .in. status), &
        "std: 2d input of shape [2, 1]")

    s1 = -1.0
    x12(1, :) = [0.213, -234.123]
    status = NF_STATUS_UNDEFINED
    call std(x12, s1, m1, dim=1, status=status)
    call tc%assert_true (all(x12(1, :) == m1) .and. all(s1 == 0.0_PREC) .and. (NF_STATUS_OK .in. status), &
        "std: 2d input of shape [1, 2]")

end subroutine



subroutine test_mean_std_input_checks (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (status_t) :: status
    real (PREC), dimension(1) :: x
    real (PREC), dimension(:,:), allocatable :: x2d
    real (PREC), dimension(10) :: mean_x, std_x
    integer :: m, n

    tc => tests%add_test ('Input checks for MEAN, STD')

    ! --- 2d API: Invalid DIM ---

    m = 2
    n = 2
    allocate (x2d(m,n))

    status = NF_STATUS_UNDEFINED
    call std (x2d, s=std_x(1:m), dim=3, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 2d API: Invalid DIM > 2')

    status = NF_STATUS_UNDEFINED
    call std (x2d, s=std_x(1:m), dim=0, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 2d API: Invalid DIM < 1')

    deallocate (x2d)

    ! --- 2d API: Non-conformable MEAN_X, STD_x arguments ---

    ! Mean
    m = 4
    n = 2
    allocate (x2d(m, n))

    status = NF_STATUS_UNDEFINED
    call mean (x2d, m=mean_x(1:n-1), dim=1, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'MEAN, 2d API: Invalid size(MEAN_X), dim=1')

    status = NF_STATUS_UNDEFINED
    call mean (x2d, m=mean_x(1:m+1), dim=2, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'MEAN, 2d API: Invalid size(MEAN_X), dim=2')

    ! Std
    status = NF_STATUS_UNDEFINED
    call std (x2d, s=std_x(1:n+1), dim=1, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 2d API: Invalid size(STD_X), dim=1')

    status = NF_STATUS_UNDEFINED
    call std (x2d, s=std_x(1:m-1), dim=2, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 2d API: Invalid size(STD_X), dim=2')

    ! Std: Invalid MEAN_X
    status = NF_STATUS_UNDEFINED
    call std (x2d, s=std_x(1:n), m=mean_x(1:n+1), dim=1, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 2d API: Invalid size(MEAN_X), dim=1')

    status = NF_STATUS_UNDEFINED
    call std (x2d, s=std_x(1:m), m=mean_x(1:m+1), dim=2, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 2d API: Invalid size(MEAN_X), dim=2')

    deallocate (x2d)

    ! --- 1d API: Invalid DOF ---

    status = NF_STATUS_UNDEFINED
    call std (x, s=std_x(1), dof=-1, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 1d API: Invalid DOF < 0')

    status = NF_STATUS_UNDEFINED
    call std (x, s=std_x(1), dof=2, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 1d API: Invalid DOF > 1')


    ! --- 2d API: Invalid DOF ---

    m = 5
    n = 2
    allocate (x2d(m,n))

    status = NF_STATUS_UNDEFINED
    call std (x2d, s=std_x(1:m), dim=1, dof=-1, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 2d API: Invalid DOF < 0')

    status = NF_STATUS_UNDEFINED
    call std (x2d, s=std_x(1:m), dim=1, dof=3, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'STD, 2d API: Invalid DOF > 1')

    deallocate (x2d)

end subroutine



subroutine test_mean_std_1d (tests)

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



subroutine test_mean_std_2d (tests)
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



subroutine test_standardize_1d (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(10) :: x, x_ok, x_orig
    real (PREC), dimension(0) :: x0
    real (PREC) :: std_x, mean_x, shift_x, scale_x, mean_ok, std_ok
    logical :: values_ok
    type (status_t) :: status

    tc => tests%add_test ('STANDARDIZE 1d API')

    call set_seed (12345)

    call random_number (x_orig)
    x = x_orig
    x_ok = x

    call std (x, m=mean_ok, s=std_ok)

    ! --- Degenerate input ---

    call standardize (x0, mean_x=mean_x, std_x=std_x, status=status)

    ! --- No operation ---

    status = NF_STATUS_UNDEFINED
    call standardize (x, center=.false., scale=.false., mean_x=mean_x, &
        std_x=std_x, shift_x=shift_x, scale_x=scale_x, status=status)
    values_ok = all(x == x_orig) .and. mean_x==mean_ok .and. std_x==std_ok &
        .and. shift_x==0.0_PREC .and. scale_x==1.0_PREC
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, CENTER=.FALSE., SCALE=.FALSE.')

    status = NF_STATUS_UNDEFINED
    call standardize (x, center=.false., scale=.false., status=status)
    values_ok = all(x == x_orig)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, CENTER=.FALSE., SCALE=.FALSE.')

    ! --- Center ---

    x = x_orig
    x_ok = x - mean_ok
    status = NF_STATUS_UNDEFINED
    call standardize (x, center=.true., scale=.false., mean_x=mean_x, &
        std_x=std_x, shift_x=shift_x, scale_x=scale_x, status=status)
    values_ok = all(x == x_ok) .and. mean_x==mean_ok .and. std_x==std_ok &
        .and. shift_x==mean_x .and. scale_x==1.0_PREC
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, CENTER=.TRUE., SCALE=.FALSE.')

    status = NF_STATUS_UNDEFINED
    x = x_orig
    x_ok = x - mean_ok
    call standardize (x, center=.true., scale=.false., status=status)
    values_ok = all(x == x_ok)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, CENTER=.TRUE., SCALE=.FALSE.')

    ! --- Scale ---

    x = x_orig
    x_ok = x / std_ok
    status = NF_STATUS_UNDEFINED
    call standardize (x, center=.false., scale=.true., mean_x=mean_x, &
        std_x=std_x, shift_x=shift_x, scale_x=scale_x, status=status)
    values_ok = all(x == x_ok) .and. mean_x==mean_ok .and. std_x==std_ok &
        .and. shift_x==0.0_PREC .and. scale_x==std_ok
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, CENTER=.FALSE., SCALE=.TRUE.')

    x = x_orig
    x_ok = x / std_ok
    status = NF_STATUS_UNDEFINED
    call standardize (x, center=.false., scale=.true., status=status)
    values_ok = all(x == x_ok)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, CENTER=.FALSE., SCALE=.TRUE.')

    ! --- Center and scale ---

    x = x_orig
    x_ok = (x - mean_ok) / std_ok
    status = NF_STATUS_UNDEFINED
    call standardize (x, center=.true., scale=.true., mean_x=mean_x, &
        std_x=std_x, shift_x=shift_x, scale_x=scale_x, status=status)
    values_ok = all(x == x_ok) .and. mean_x==mean_ok .and. std_x==std_ok &
        .and. shift_x==mean_ok .and. scale_x==std_ok
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, CENTER=.TRUE., SCALE=.TRUE.')

    x = x_orig
    x_ok = (x - mean_ok) / std_ok
    status = NF_STATUS_UNDEFINED
    call standardize (x, center=.true., scale=.true., status=status)
    values_ok = all(x == x_ok)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, CENTER=.TRUE., SCALE=.TRUE.')

end subroutine



subroutine test_standardize_2d (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    integer, parameter :: m = 7, n = 4
    real (PREC), dimension(m,n) :: x, x_ok, x_orig
    real (PREC), dimension(max(m,n)) :: std_x, mean_x, shift_x, scale_x, mean_ok, std_ok
    integer :: Nvars, Nobs
    logical :: values_ok
    integer :: k, dim, i
    type (status_t) :: status

    tc => tests%add_test ('STANDARDIZE 2d API')

    call set_seed (12345)

    call random_number (x_orig)
    x = x_orig
    x_ok = x

    ! --- Degenerate inputs ---

    ! Zero obs, positive number of variables
    ! This should be considered as an invalid argument since the
    ! mean and std. dev. are not defined for 0 obs.
    Nvars = min(m,n)
    Nobs = 0
    dim = 1
    status = NF_STATUS_UNDEFINED
    call standardize (x(1:Nobs,1:Nvars), dim=dim, mean_x=mean_x(1:Nvars), &
        std_x=std_x(1:Nvars), status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Nobs = 0 , Nvar > 0, dim = 1')

    dim = 2
    status = NF_STATUS_UNDEFINED
    call standardize (x(1:Nvars,1:Nobs), dim=dim, mean_x=mean_x(1:Nvars), &
        std_x=std_x(1:Nvars), status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Nobs = 0 , Nvar > 0, dim = 2')

    ! Zero vars, positive number of variables.
    ! This should return NF_STATUS_OK.
    Nvars = 0
    Nobs = min(m,n)
    dim = 1
    status = NF_STATUS_UNDEFINED
    call standardize (x(1:Nobs,1:Nvars), dim=dim, mean_x=mean_x(1:Nvars), &
        std_x=std_x(1:Nvars), status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'Nobs > 0 , Nvar = 0, dim = 1')

    dim = 2
    status = NF_STATUS_UNDEFINED
    call standardize (x(1:Nvars,1:Nobs), dim=dim, mean_x=mean_x(1:Nvars), &
        std_x=std_x(1:Nvars), status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'Nobs > 0 , Nvar = 0, dim = 2')

    ! --- No operation ---

    ! DIM = 1
    dim = 1
    k = n
    call std (x, dim=dim, m=mean_ok(1:k), s=std_ok(1:k))

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.false., scale=.false., mean_x=mean_x(1:k), &
        std_x=std_x(1:k), shift_x=shift_x(1:k), scale_x=scale_x(1:k), status=status)
    values_ok = all(x == x_orig) .and. all(mean_x(1:k) == mean_ok(1:k)) &
        .and. all(std_x(1:k) == std_ok(1:k)) &
        .and. all(shift_x(1:k) == 0.0_PREC) &
        .and. all(scale_x(1:k) == 1.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, DIM=1, CENTER=.FALSE., SCALE=.FALSE.')

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.false., scale=.false., status=status)
    values_ok = all(x == x_orig)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, DIM=1, CENTER=.FALSE., SCALE=.FALSE.')

    ! DIM = 2
    k = m
    dim = 2
    call std (x, dim=dim, m=mean_ok(1:k), s=std_ok(1:k))

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.false., scale=.false., &
        mean_x=mean_x(1:k), std_x=std_x(1:k), shift_x=shift_x(1:k), &
        scale_x=scale_x(1:k), status=status)

    values_ok = all(x == x_orig) .and. all(mean_x(1:k) == mean_ok(1:k)) &
        .and. all(std_x(1:k) == std_ok(1:k)) &
        .and. all(shift_x(1:k) == 0.0_PREC) &
        .and. all(scale_x(1:k) == 1.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, DIM=1, CENTER=.FALSE., SCALE=.FALSE.')

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.false., scale=.false., status=status)

    values_ok = all(x == x_orig)

    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, DIM=1, CENTER=.FALSE., SCALE=.FALSE.')

    ! --- Center ---

    ! DIM = 1
    dim = 1
    k = n
    mean_ok = 0.0
    std_ok = 0.0
    call std (x, dim=dim, m=mean_ok(1:k), s=std_ok(1:k))
    x = x_orig
    do i = 1, n
        x_ok(:,i) = x_orig(:,i) - mean_ok(i)
    end do

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.true., scale=.false., &
        mean_x=mean_x(1:k), std_x=std_x(1:k), shift_x=shift_x(1:k), &
        scale_x=scale_x(1:k), status=status)

    values_ok = all(x == x_ok) .and. all(mean_x(1:k) == mean_ok(1:k)) &
        .and. all(std_x(1:k) == std_ok(1:k)) &
        .and. all(shift_x(1:k) == mean_ok(1:k)) &
        .and. all(scale_x(1:k) == 1.0_PREC)

    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, DIM=1, CENTER=.TRUE., SCALE=.FALSE.')

    x = x_orig

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.true., scale=.false., status=status)
    values_ok = all(x == x_ok)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, DIM=1, CENTER=.TRUE., SCALE=.FALSE.')

    ! DIM = 2
    x = x_orig
    dim = 2
    k = m
    mean_ok = 0.0
    std_ok = 0.0
    call std (x, dim=dim, m=mean_ok(1:k), s=std_ok(1:k))
    do i = 1, m
        x_ok(i,:) = x_orig(i,:) - mean_ok(i)
    end do

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.true., scale=.false., &
        mean_x=mean_x(1:k), std_x=std_x(1:k), shift_x=shift_x(1:k), &
        scale_x=scale_x(1:k), status=status)

    values_ok = all(x == x_ok) .and. all(mean_x(1:k) == mean_ok(1:k)) &
        .and. all(std_x(1:k) == std_ok(1:k)) &
        .and. all(shift_x(1:k) == mean_ok(1:k)) &
        .and. all(scale_x(1:k) == 1.0_PREC)

    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, DIM=2, CENTER=.TRUE., SCALE=.FALSE.')

    x = x_orig

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.true., scale=.false., status=status)

    values_ok = all(x == x_ok)

    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, DIM=2, CENTER=.TRUE., SCALE=.FALSE.')

    ! --- Scale ---

    ! DIM = 1
    x = x_orig
    dim = 1
    k = n
    mean_ok = 0.0
    std_ok = 0.0
    call std (x, dim=dim, m=mean_ok(1:k), s=std_ok(1:k))
    do i = 1, n
        x_ok(:,i) = x_orig(:,i) / std_ok(i)
    end do

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.false., scale=.true., &
        mean_x=mean_x(1:k), std_x=std_x(1:k), shift_x=shift_x(1:k), &
        scale_x=scale_x(1:k), status=status)

    values_ok = all(x == x_ok) .and. all(mean_x(1:k) == mean_ok(1:k)) &
        .and. all(std_x(1:k) == std_ok(1:k)) &
        .and. all(shift_x(1:k) == 0.0_PREC) &
        .and. all(scale_x(1:k) == std_ok(1:k))

    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, DIM=1, CENTER=.FALSE., SCALE=.TRUE.')

    x = x_orig

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.false., scale=.true., status=status)
    values_ok = all(x == x_ok)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, DIM=1, CENTER=.FALSE., SCALE=.TRUE.')

    ! DIM = 2
    x = x_orig
    dim = 2
    k = m
    mean_ok = 0.0
    std_ok = 0.0
    call std (x, dim=dim, m=mean_ok(1:k), s=std_ok(1:k))
    do i = 1, m
        x_ok(i,:) = x_orig(i,:) / std_ok(i)
    end do

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.false., scale=.true., &
        mean_x=mean_x(1:k), std_x=std_x(1:k), shift_x=shift_x(1:k), &
        scale_x=scale_x(1:k), status=status)

    values_ok = all(x == x_ok) .and. all(mean_x(1:k) == mean_ok(1:k)) &
        .and. all(std_x(1:k) == std_ok(1:k)) &
        .and. all(shift_x(1:k) == 0.0_PREC) &
        .and. all(scale_x(1:k) == std_ok(1:k))

    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, DIM=1, CENTER=.FALSE., SCALE=.TRUE.')

    x = x_orig

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.false., scale=.true., status=status)
    values_ok = all(x == x_ok)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, DIM=1, CENTER=.FALSE., SCALE=.TRUE.')

    ! --- Center and scale ---

    ! DIM = 1
    x = x_orig
    dim = 1
    k = n
    mean_ok = 0.0
    std_ok = 0.0
    call std (x, dim=dim, m=mean_ok(1:k), s=std_ok(1:k))
    do i = 1, n
        x_ok(:,i) = (x_orig(:,i) - mean_ok(i)) / std_ok(i)
    end do

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.true., scale=.true., &
        mean_x=mean_x(1:k), std_x=std_x(1:k), shift_x=shift_x(1:k), &
        scale_x=scale_x(1:k), status=status)

    values_ok = all(x == x_ok) .and. all(mean_x(1:k) == mean_ok(1:k)) &
        .and. all(std_x(1:k) == std_ok(1:k)) &
        .and. all(shift_x(1:k) == mean_ok(1:k)) &
        .and. all(scale_x(1:k) == std_ok(1:k))

    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, DIM=1, CENTER=.TRUE., SCALE=.TRUE.')

    x = x_orig

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.true., scale=.true., status=status)
    values_ok = all(x == x_ok)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, DIM=1, CENTER=.TRUE., SCALE=.TRUE.')

    ! DIM = 2
    x = x_orig
    dim = 2
    k = m
    mean_ok = 0.0
    std_ok = 0.0
    call std (x, dim=dim, m=mean_ok(1:k), s=std_ok(1:k))
    do i = 1, m
        x_ok(i,:) = (x_orig(i,:) - mean_ok(i)) / std_ok(i)
    end do

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.true., scale=.true., &
        mean_x=mean_x(1:k), std_x=std_x(1:k), shift_x=shift_x(1:k), &
        scale_x=scale_x(1:k), status=status)

    values_ok = all(x == x_ok) .and. all(mean_x(1:k) == mean_ok(1:k)) &
        .and. all(std_x(1:k) == std_ok(1:k)) &
        .and. all(shift_x(1:k) == mean_ok(1:k)) &
        .and. all(scale_x(1:k) == std_ok(1:k))

    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'All args present, DIM=1, CENTER=.TRUE., SCALE=.TRUE.')

    x = x_orig

    status = NF_STATUS_UNDEFINED
    call standardize (x, dim=dim, center=.true., scale=.true., status=status)
    values_ok = all(x == x_ok)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'No optional args present, DIM=1, CENTER=.FALSE., SCALE=.TRUE.')

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
    nobs = 100000
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
