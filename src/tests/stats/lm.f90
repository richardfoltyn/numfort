

program test_stats_lm
    !*  Unit tests for OLS estimation routines in numfort_stats_lm module.

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_testing
    use numfort_stats, only: set_seed
    use numfort_stats_lm, lm_result => lm_result_real64, lm_config => lm_config_real64
    use numfort_stats_data_helpers
    use numfort_stats_dnorm

    use blas95, only: BLAS95_GEMM => GEMM, BLAS95_GEMV => GEMV

    use fcore_testing
    use fcore_strings

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()


    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("numfort_stats_lm unit tests")

    call test_OLS_1d_args (tests)
    call test_OLS_2d_args (tests)
    call test_OLS_1d (tests)
    call test_OLS_2d (tests)
    call test_OLS_fit_equiv (tests)

    call test_predict_equiv_1d (tests)
    call test_predict_equiv_2d (tests)

    call test_pcr_equiv_1d (tests)
    call test_pcr_equiv_2d (tests)

    call test_pcr_ols_equiv_1d (tests)
    call test_pcr_ols_equiv_2d (tests)

    call tests%print ()

end subroutine



subroutine test_OLS_1d_args (tests)
    !*  Unit tests for OLS argument checking
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: x
    real (PREC), dimension(:), allocatable :: y
    real (PREC), dimension(:), allocatable :: beta

    integer :: nrhs, nobs, ncoefs
    type (lm_config) :: conf
    type (status_t) :: status

    tc => tests%add_test ("OLS[1d] argument checking")

    ! === incompatible X, Y arguments ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs+1), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, Y")

    deallocate (y, x, beta)

    ! === incompatible X^T, Y ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs, nobs), y(nobs+1), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(trans_x=.true.)
    call ols (conf, x, y, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, Y")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs-2))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===

    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs+2,nobs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(trans_x=.true.)
    call ols (conf, x, y, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS with ADD_INTERCEPT=.FALSE. ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(add_intercept=.false.)
    call ols (conf, x, y, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS with ADD_INTERCEPT=.FALSE.")

    deallocate (x, y, beta)

    ! === NOBS < NCOEFS ===
    nrhs = 3
    nobs = nrhs
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(add_intercept=.true.)
    call ols (conf, x, y, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "NOBS < NCOEFS")

    deallocate (x, y, beta)

    ! === Invalid RCOND ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(rcond=-0.1_PREC)
    call ols (conf, x, y, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Invalid RCOND")

    deallocate (x, y, beta)

end subroutine


subroutine test_OLS_2d_args (tests)
    !*  Unit tests for OLS argument checking
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: x, y
    real (PREC), dimension(:,:), allocatable :: beta

    integer :: nrhs, nlhs, nobs, ncoefs
    type (lm_config) :: conf
    type (status_t) :: status

    tc => tests%add_test ("OLS[2d] argument checking")

    ! === incompatible X, Y arguments ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs+1,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, Y")

    deallocate (x, y, beta)

    ! === incompatible X^T, Y ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs, nobs), y(nobs+1,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x, y, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, Y")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(nrhs-1,nlhs))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs+2,nobs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x, y, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS with ADD_INTERCEPT=.FALSE. ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS with ADD_INTERCEPT=.FALSE.")

    deallocate (x, y, beta)

    ! === incompabible Y, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs-1))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays Y, COEFS")

    deallocate (x, y, beta)

    ! === NOBS < NCOEFS ===
    nrhs = 3
    nlhs = 2
    nobs = nrhs
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "NOBS < NCOEFS")

    deallocate (x, y, beta)

    ! === Invalid RCOND ===
    nrhs = 3
    nlhs = 2
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config (rcond=-1.0_PREC)
    call ols (conf, x, y, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Invalid RCOND")

    deallocate (x, y, beta)

end subroutine



subroutine test_OLS_2d (tests)
    !*  Unit tests for various OLS regression problems using the 2d API
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: x, y, beta, beta_ok, work
    integer :: nrhs, nlhs, ncoefs, nobs
    integer :: i
    integer :: rank
    real (PREC), parameter :: RTOL = 1.0e-4_PREC, ATOL=1.0e-12_PREC
        ! Arguments for ALL_CLOSE()
    real (PREC) :: rsq
    logical :: all_ok, values_ok
    type (lm_config) :: conf
    type (lm_result), allocatable :: lm
    type (status_t) :: status

    tc => tests%add_test ('OLS[2d] unit tests')

    call set_seed (123)

    ! === exactly determined system with NOBS = NRHS, no const ===
    nrhs = 10
    nobs = 10
    ncoefs = nrhs
    nlhs = 3

    allocate (x(nobs,nrhs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (x)
    call random_number (beta_ok)

    call BLAS95_GEMM (x, beta_ok, y)

    status = NF_STATUS_UNDEFINED
    allocate (lm)
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, rank=rank, res=lm, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        'Exactly identified system: estimation')
    deallocate (lm)

    ! Test without optional COEFS argument
    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.false.)
    allocate (lm)
    call ols (conf, x, y, rank=rank, res=lm, status=status)
    values_ok = all_close (lm%coefs_multi, beta_ok, rtol=RTOL, atol=ATOL)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Exactly identified system: estimation w/o COEFS argument')

    deallocate (x, y, beta, beta_ok)
    deallocate (lm)

    ! === overdetermined system without intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs
    nlhs = 2

    allocate (x(nobs,nrhs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    call BLAS95_GEMM (x, beta_ok, y)

    status = NF_STATUS_OK
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        "Underdetermined system")

    deallocate (x, y, beta, beta_ok)


    ! === overdetermined system with intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs + 1
    nlhs = 11

    allocate (x(nobs,nrhs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    forall (i=1:nlhs) y(:,i) = beta_ok(1,i)
    allocate (work(nrhs,nlhs), source=beta_ok(2:,:))
    call BLAS95_GEMM (x, work, y, beta=1.0_PREC)
    deallocate (work)

    status = NF_STATUS_OK
    conf = lm_config (add_intercept=.true.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    ! Apply less strict check on beta, value of 1.0e-3 does not seem to work
    ! for NLHS=11
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=1.0e-2_PREC, atol=ATOL), &
        "Underdetermined system, ADD_INTERCEPT=.TRUE.")

    deallocate (x, y, beta, beta_ok)


    ! === overdetermined system with intercept, transposed X ===
    nrhs = 3
    nobs = 500
    ncoefs = nrhs + 1
    nlhs = 7

    allocate (x(nrhs,nobs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    forall (i=1:nlhs) y(:,i) = beta_ok(1,i)
    allocate (work(nrhs,nlhs), source=beta_ok(2:,:))
    call BLAS95_GEMM (x, work, y, beta=1.0_PREC, transa='T')
    deallocate (work)

    status = NF_STATUS_OK
    conf = lm_config (add_intercept=.true., trans_x=.true.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=1.0e-2_PREC, atol=ATOL), &
        "Underdetermined system, ADD_INTERCEPT=.TRUE., TRANS_X=.TRUE.")

    deallocate (x, y, beta, beta_ok)

end subroutine



subroutine test_OLS_1d (tests)
    !*  Unit tests for various OLS regression problems using the 1d API
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: x
    real (PREC), dimension(:), allocatable :: y, beta, beta_ok, work
    integer :: nrhs, ncoefs, nobs
    integer :: rank
    real (PREC), parameter :: RTOL = 1.0e-4_PREC, ATOL=1.0e-12_PREC
        ! Arguments for ALL_CLOSE()
    real (PREC) :: rsq
    type (status_t) :: status
    type (lm_config) :: conf
    type (lm_result), allocatable :: lm

    tc => tests%add_test ('OLS[1d] unit tests')

    call set_seed (456)

    ! === exactly determined system with NOBS = NRHS, no const ===

    nrhs = 10
    nobs = 10
    ncoefs = nrhs

    call random_sample (nrhs, nobs, x, y, coefs=beta_ok, add_intercept=.false.)

    status = NF_STATUS_UNDEFINED
    allocate (lm)
    allocate (beta(ncoefs))
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, rank=rank, res=lm, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        'Exactly identified system: estimation')

    ! Test w/o optional COEFS argument, coefs are stored in LM_DATA
    deallocate (lm)
    allocate (lm)
    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, rank=rank, res=lm, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (lm%coefs, beta_ok, rtol=RTOL, atol=ATOL), &
        'Exactly identified system: estimation w/o COEFS argument')

    ! Compute R^2 which must be 1.0 in this case
    status = NF_STATUS_UNDEFINED
    call post_estim (lm, x, y, rsq=rsq, status=status)
    call tc%assert_true (abs(rsq-1.0_PREC) < 1.0e-12 .and. status == NF_STATUS_OK, &
        'Exactly identified system: R^2')

    deallocate (x, y, beta, beta_ok)

    ! === overdetermined system without intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_ok(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    call BLAS95_GEMV (x, beta_ok, y)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        "Overdetermined system")

    deallocate (x, y, beta, beta_ok)


    ! === overdetermined system with intercept ===

    nrhs = 5
    nobs = 500
    ncoefs = nrhs + 1

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_ok(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    y(:) = beta_ok(1)
    allocate (work(nrhs), source=beta_ok(2:))
    call BLAS95_GEMV (x, work, y, beta=1.0_PREC)
    deallocate (work)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.true.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        "Overdetermined system, ADD_INTERCEPT=.TRUE.")

    deallocate (x, y, beta, beta_ok)


    ! === overdetermined system with intercept, transposed X ===
    nrhs = 2
    nobs = 500
    ncoefs = nrhs + 1

    allocate (x(nrhs,nobs), y(nobs), beta(ncoefs), beta_ok(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    y(:) = beta_ok(1)
    allocate (work(nrhs), source=beta_ok(2:))
    call BLAS95_GEMV (x, work, y, beta=1.0_PREC, trans='T')
    deallocate (work)

    status = NF_STATUS_OK
    conf = lm_config (add_intercept=.true., trans_x=.true.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        "Overdetermined system, ADD_INTERCEPT=.TRUE., TRANS_X=.TRUE.")

    deallocate (x, y, beta, beta_ok)

end subroutine



subroutine test_OLS_fit_equiv (tests)
    !*  Unit tests for equivalent formulations of the same LS problem,
    !   with OLS() being called with non-transposed and transposed
    !   input and output data.
    class (test_suite) :: tests

    type (test_case), pointer :: tc
    type (lm_config) :: conf
    type (lm_result), allocatable :: res
    real (PREC), dimension(:,:), allocatable :: X, X_T, Y, Y_T
    real (PREC), dimension(:,:), allocatable :: coefs, coefs_T, coefs_true
    real (PREC), dimension(:), allocatable ::  intercept, intercept_true
    integer :: Nobs, Nrhs, Nconst, Ncoefs, Nlhs
    type (status_t) :: status
    logical :: values_ok
    real (PREC), parameter :: atol = 1.0e-12_PREC, rtol = 1.0e-4_PREC
    integer :: jlhs

    tc => tests%add_test ('OLS fitting: equivalent problems')

    call set_seed (1234)

    ! Test data set

    Nobs = 11
    Nrhs = 3
    Nlhs = 5

    allocate (X(Nobs,Nrhs))
    allocate (coefs_true(Nrhs,Nlhs))
    allocate (intercept_true(Nlhs))

    call random_number (X)
    call random_number (coefs_true)
    call random_number (intercept_true)

    coefs_true(:,:) = (coefs_true - 0.5d0) * 5.0d0
    intercept_true(:) = (intercept_true - 0.5d0) * 10.d0

    allocate (Y(Nobs,Nlhs))

    Y(:,:) = spread (intercept_true, dim=1, ncopies=Nobs)

    call BLAS95_GEMM (X, coefs_true, Y, alpha=1.0_PREC, beta=1.0_PREC)

    allocate (X_T(Nrhs,Nobs), Y_T(Nlhs,Nobs))
    X_T(:,:) = transpose(X)
    Y_T(:,:) = transpose(Y)

    ! --- 1d API, no intercept in COEFS ---

    ! trans_x = .false.
    allocate (coefs(Nrhs,1), source=huge(0.0_PREC))
    allocate (intercept(1), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    call ols (x=X, y=Y(:,1), coefs=coefs(:,1), intercept=intercept(1), res=res, status=status)

    values_ok = all_close (coefs(:,1), coefs_true(:,1), atol=atol, rtol=rtol) &
        .and. all_close (intercept(1), intercept_true(1), atol=atol, rtol=rtol) &
        .and. all_close (res%coefs, coefs_true(:,1), atol=atol, rtol=rtol) &
        .and. all_close (res%intercept, intercept_true(1), atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_X=.FALSE., no intercept in coefs')

    deallocate (res)

    ! trans_x = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x=X_T, y=Y(:,1), coefs=coefs(:,1), intercept=intercept(1), &
        res=res, status=status)

    values_ok = all_close (coefs(:,1), coefs_true(:,1), atol=atol, rtol=rtol) &
        .and. all_close (intercept(1), intercept_true(1), atol=atol, rtol=rtol) &
        .and. all_close (res%coefs, coefs_true(:,1), atol=atol, rtol=rtol) &
        .and. all_close (res%intercept, intercept_true(1), atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_X=.TRUE., no intercept in coefs')

    deallocate (res)
    deallocate (coefs, intercept)

    ! --- 1d API, add intercept to COEFS ---

    ! trans_x = .false.
    Nconst = 1
    Ncoefs = Nrhs + Nconst
    allocate (coefs(Ncoefs,1), source=huge(0.0_PREC))
    allocate (intercept(1), source=huge(0.0_PREC))
    allocate (res)

    jlhs = 2
    status = NF_STATUS_UNDEFINED
    call ols (x=X, y=Y(:,jlhs), coefs=coefs(:,1), intercept=intercept(1), res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,1), coefs_true(:,jlhs), atol=atol, rtol=rtol) &
        .and. all_close (intercept(1), intercept_true(jlhs), atol=atol, rtol=rtol) &
        .and. all_close (res%coefs, coefs_true(:,jlhs), atol=atol, rtol=rtol) &
        .and. all_close (res%intercept, intercept_true(jlhs), atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_X=.FALSE., intercept in coefs')

    deallocate (res)

    ! trans_x = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x=X_T, y=Y(:,jlhs), coefs=coefs(:,1), intercept=intercept(1), &
        res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,1), coefs_true(:,jlhs), atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,1), intercept_true(jlhs), atol=atol, rtol=rtol) &
        .and. all_close (intercept(1), intercept_true(jlhs), atol=atol, rtol=rtol) &
        .and. all_close (res%coefs, coefs_true(:,jlhs), atol=atol, rtol=rtol) &
        .and. all_close (res%intercept, intercept_true(jlhs), atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_X=.TRUE., intercept in coefs')

    deallocate (res)
    deallocate (coefs, intercept)

    ! --- 2d API, no intercept in COEFS ---

    ! trans_x = .false.
    allocate (coefs(Nrhs,Nlhs), source=huge(0.0_PREC))
    allocate (intercept(Nlhs), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    call ols (x=X, y=Y, coefs=coefs, intercept=intercept, res=res, status=status)

    values_ok = all_close (coefs, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '2d API: TRANS_X=.FALSE., no intercept in coefs')

    deallocate (res)

    ! trans_x = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x=X_T, y=Y, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_X=.TRUE., no intercept in coefs')

    deallocate (res)

    ! trans_x = .false., trans_y = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.false., trans_y=.true.)
    call ols (conf, x=X, y=Y_T, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_Y=.TRUE., no intercept in coefs')

    deallocate (res)

    ! trans_x = .true., trans_y = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true., trans_y=.true.)
    call ols (conf, x=X_T, y=Y_T, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_X=.TRUE., TRANS_Y=.TRUE., no intercept in coefs')

    deallocate (res)
    deallocate (coefs, intercept)


    ! trans_coefs = .true.
    allocate (coefs_T(Nlhs,Nrhs), source=huge(0.0_PREC))
    allocate (intercept(Nlhs), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config(trans_coefs=.true.)
    call ols (conf, x=X, y=Y, coefs=coefs_T, intercept=intercept, res=res, status=status)

    values_ok = all_close (transpose(coefs_T), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '2d API: TRANS_COEFS=.TRUE., no intercept in coefs')

    deallocate (coefs_T, intercept)
    deallocate (res)

    ! --- 2d API, add intercept to COEFS ---

    ! trans_x = .false.
    Nconst = 1
    Ncoefs = Nrhs + Nconst
    allocate (coefs(Ncoefs,Nlhs), source=huge(0.0_PREC))
    allocate (intercept(Nlhs), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    call ols (x=X, y=Y, coefs=coefs, intercept=intercept, res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,:), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,:), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_X=.FALSE., intercept in coefs')

    deallocate (res)

    ! trans_x = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x=X_T, y=Y, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,:), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,:), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_X=.TRUE., intercept in coefs')

    deallocate (res)

    ! trans_y = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.false., trans_y=.true.)
    call ols (conf, x=X, y=Y_T, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,:), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,:), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_Y=.TRUE., intercept in coefs')

    deallocate (res)


    ! trans_x = .true., trans_y = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true., trans_y=.true.)
    call ols (conf, x=X_T, y=Y_T, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,:), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,:), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_X=.TRUE., TRANS_Y=.TRUE., intercept in coefs')

    deallocate (res)
    deallocate (coefs, intercept)


    ! trans_coefs = .true.
    Nconst = 1
    Ncoefs = Nrhs + Nconst
    allocate (coefs_T(Nlhs,Ncoefs), source=huge(0.0_PREC))
    allocate (intercept(Nlhs), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_coefs=.true.)
    call ols (conf, x=X, y=Y, coefs=coefs_T, intercept=intercept, res=res, status=status)

    values_ok = all_close (transpose(coefs_T(:,1+Nconst:Ncoefs)), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs_T(:,1), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        '1d API: TRANS_COEFS=.TRUE., intercept in coefs')

    deallocate (coefs_T, intercept)
    deallocate (res)

end subroutine



subroutine test_predict_equiv_1d (tests)
    !*  Test equivalent calls to PREDICT that use transposed variations of
    !   input and output data
    class (test_suite) :: tests

    type (test_case), pointer :: tc
    type (lm_config) :: conf
    type (lm_result) :: res
    real (PREC), dimension(:,:), allocatable :: X, Xpred, Xpred_T
    real (PREC), dimension(:), allocatable :: Y, Ypred
    real (PREC), dimension(:), allocatable :: Ypred_true
    real (PREC), dimension(:), allocatable :: coefs
    real (PREC), dimension(:), allocatable :: coefs_intercept
    real (PREC) ::  intercept
    integer :: Nobs, Nrhs, Nconst, Ncoefs, Nlhs, Nobs_pred
    type (status_t) :: status
    logical :: values_ok
    real (PREC), parameter :: atol = 1.0e-12_PREC, rtol = 1.0e-4_PREC
    integer :: jlhs
    logical :: trans_x, trans_y, trans_coefs

    tc => tests%add_test ('PREDICT[1d]: equivalent problems')

    call set_seed (45678)

    ! --- Training data set ---

    Nobs = 11
    Nrhs = 3
    Nconst = 1
    Ncoefs = Nrhs + Nconst

    allocate (X(Nobs,Nrhs))
    allocate (coefs(Nrhs))

    call random_number (X)
    call random_number (coefs)
    call random_number (intercept)

    coefs(:) = (coefs - 0.5d0) * 5.0d0
    intercept = (intercept - 0.5d0) * 10.d0

    allocate (Y(Nobs))

    Y(:) = intercept

    call BLAS95_GEMV (X, coefs, Y, alpha=1.0_PREC, beta=1.0_PREC)

    ! Coefs + intercept in one array
    allocate (coefs_intercept(Ncoefs))
    coefs_intercept(1) = intercept
    coefs_intercept(1+Nconst:Ncoefs) = coefs

    ! Create OLS result objects so we can test PREDICT(LM_RESULT,...)
    call ols (X=X, Y=Y, res=res, status=status)

    ! --- Test data ---
    Nobs_pred = Nobs + 12

    allocate (Xpred(Nobs_pred,Nrhs), Xpred_T(Nrhs,Nobs_pred))
    allocate (Ypred_true(Nobs_pred))

    call random_number (Xpred)
    Xpred_T(:,:) = transpose (Xpred)

    Ypred_true(:) = intercept
    call BLAS95_GEMV (Xpred, coefs, Ypred_true, alpha=1.0_PREC, beta=1.0_PREC)

    ! Arrays where PREDICT writes its results
    allocate (Ypred(Nobs_pred))

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred, coefs, Ypred, intercept, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'TRANS_X=.FALSE.')

    ! trans_x = .false., called via PREDICT(LM_RESULT,...)
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    Ypred(:) = huge(0.0_PREC)
    call predict (res, Xpred, Ypred, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'TRANS_X=.FALSE., using LM_RESULT')

    ! trans_x = .false., extended COEFS, INTERCEPT present
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred, &
        intercept=intercept, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'TRANS_X=.FALSE., extended COEFS, INTERCEPT present')

    ! trans_x = .false., extended COEFS, INTERCEPT not present
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'TRANS_X=.FALSE., extended COEFS')

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x = .true.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs, Ypred, intercept, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'TRANS_X=.FALSE.')

    ! trans_x = .true., using LM_RESULT
    status = NF_STATUS_UNDEFINED
    trans_x = .true.
    Ypred(:) = huge(0.0_PREC)
    call predict (res, Xpred_T, Ypred, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'TRANS_X=.FALSE., using LM_RESULT')

    ! trans_x = .true., extended COEFS, INTERCEPT present
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred, &
        intercept=intercept, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'TRANS_X=.TRUE., extended COEFS, INTERCEPT present')

    ! trans_x = .true., extended COEFS, INTERCEPT not present
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'TRANS_X=.TRUE., extended COEFS')

end subroutine



subroutine test_predict_equiv_2d (tests)
    !*  Test equivalent calls to PREDICT that use transposed variations of
    !   input and output data
    class (test_suite) :: tests

    type (test_case), pointer :: tc
    type (lm_config) :: conf
    type (lm_result) :: res
    real (PREC), dimension(:,:), allocatable :: X, Xpred, Xpred_T
    real (PREC), dimension(:,:), allocatable :: Y, Ypred, Ypred_T
    real (PREC), dimension(:,:), allocatable :: Ypred_true, Ypred_true_T
    real (PREC), dimension(:,:), allocatable :: coefs, coefs_T
    real (PREC), dimension(:,:), allocatable :: coefs_intercept, coefs_intercept_T
    real (PREC), dimension(:), allocatable ::  intercept
    integer :: Nobs, Nrhs, Nconst, Ncoefs, Nlhs, Nobs_pred
    type (status_t) :: status
    logical :: values_ok
    real (PREC), parameter :: atol = 1.0e-12_PREC, rtol = 1.0e-4_PREC
    logical :: trans_x, trans_y, trans_coefs

    tc => tests%add_test ('PREDICT[2d]: equivalent problems')

    call set_seed (7890)

    ! --- Training data set ---

    Nobs = 11
    Nrhs = 3
    Nlhs = 5
    Nconst = 1
    Ncoefs = Nrhs + Nconst

    allocate (X(Nobs,Nrhs))
    allocate (coefs(Nrhs,Nlhs))
    allocate (intercept(Nlhs))

    call random_number (X)
    call random_number (coefs)
    call random_number (intercept)

    coefs(:,:) = (coefs - 0.5d0) * 5.0d0
    intercept(:) = (intercept - 0.5d0) * 10.d0

    allocate (Y(Nobs,Nlhs))

    Y(:,:) = spread (intercept, dim=1, ncopies=Nobs)

    call BLAS95_GEMM (X, coefs, Y, alpha=1.0_PREC, beta=1.0_PREC)

    allocate (coefs_T(Nlhs,Nrhs))
    coefs_T(:,:) = transpose(coefs)

    ! Coefs + intercept in one array
    allocate (coefs_intercept(Ncoefs,Nlhs))
    coefs_intercept(1,:) = intercept
    coefs_intercept(1+Nconst:Ncoefs,:) = coefs

    allocate (coefs_intercept_T(Nlhs,Ncoefs))
    coefs_intercept_T(:,:) = transpose (coefs_intercept)

    ! Create OLS result objects so we can test PREDICT(LM_RESULT,...)
    call ols (X=X, Y=Y, res=res, status=status)

    ! --- Test data ---
    Nobs_pred = Nobs + 12

    allocate (Xpred(Nobs_pred,Nrhs), Xpred_T(Nrhs,Nobs_pred))
    allocate (Ypred_true(Nobs_pred,Nlhs), Ypred_true_T(Nlhs,Nobs_pred))

    call random_number (Xpred)
    Xpred_T(:,:) = transpose (Xpred)

    Ypred_true(:,:) = spread (intercept, dim=1, ncopies=Nobs_pred)
    call BLAS95_GEMM (Xpred, coefs, Ypred_true, alpha=1.0_PREC, beta=1.0_PREC)

    Ypred_true_T(:,:) = transpose (Ypred_true)

    ! Arrays where PREDICT writes its results
    allocate (Ypred(Nobs_pred,Nlhs), Ypred_T(Nlhs,Nobs_pred))


    ! --- Tests using LM_RESULT ---

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (res, Xpred, Ypred, trans_x=trans_x, trans_y=trans_y, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'LM_RESULT: TRANS_X=.FALSE.')

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x = .true.
    trans_y = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (res, Xpred_T, Ypred, trans_x=trans_x, trans_y=trans_y, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'LM_RESULT: TRANS_X=.TRUE.')

    ! trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (res, Xpred, Ypred_T, trans_x=trans_x, trans_y=trans_y, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'LM_RESULT: TRANS_Y=.TRUE.')

    ! trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x = .true.
    trans_y = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (res, Xpred_T, Ypred_T, trans_x=trans_x, trans_y=trans_y, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'LM_RESULT: TRANS_X=.TRUE., TRANS_Y=.TRUE.')

    ! --- Test using COEFS and INTERCEPT ---

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'COEFS: TRANS_X=.FALSE.')

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'COEFS: TRANS_X=.TRUE.')

    ! trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'COEFS: TRANS_Y=.TRUE.')

    ! trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'COEFS: TRANS_X=.TRUE., TRANS_Y=.TRUE.')

    ! trans_coefs = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_T, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'COEFS: TRANS_COEFS=.TRUE.')

    ! trans_coefs = .true., trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_T, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'COEFS: TRANS_COEFS=.TRUE., TRANS_X=.TRUE.')

    ! trans_coefs = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_T, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'COEFS: TRANS_COEFS=.TRUE., TRANS_Y=.TRUE.')

    ! trans_coefs = .true., trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_T, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'COEFS: TRANS_COEFS=.TRUE., TRANS_X=.TRUE., TRANS_Y=.TRUE.')


    ! --- Extended COEFS, INTERCEPT present ---

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Ext. COEFS, INTERCEPT: TRANS_X=.FALSE.')

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Ext. COEFS, INTERCEPT: TRANS_X=.TRUE.')

    ! trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Ext. COEFS, INTERCEPT: TRANS_Y=.TRUE.')

    ! trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Ext. COEFS, INTERCEPT: TRANS_X=.TRUE., TRANS_Y=.TRUE.')

    ! trans_coefs = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept_T, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Ext. COEFS, INTERCEPT: TRANS_COEFS=.TRUE.')

    ! trans_coefs = .true., trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept_T, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Ext. COEFS, INTERCEPT: TRANS_COEFS=.TRUE., TRANS_X=.TRUE.')

    ! trans_coefs = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept_T, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Ext. COEFS, INTERCEPT: TRANS_COEFS=.TRUE., TRANS_Y=.TRUE.')

    ! trans_coefs = .true., trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept_T, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Ext. COEFS, INTERCEPT: TRANS_COEFS=.TRUE., TRANS_X=.TRUE., TRANS_Y=.TRUE.')

    ! --- Extended COEFS, INTERCEPT not present ---

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Extended COEFS: TRANS_X=.FALSE.')

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Extended COEFS: TRANS_X=.TRUE.')

    ! trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred_T, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Extended COEFS: TRANS_Y=.TRUE.')

    ! trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred_T, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Extended COEFS: TRANS_X=.TRUE., TRANS_Y=.TRUE.')

    ! trans_coefs = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept_T, Ypred, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Extended COEFS: TRANS_COEFS=.TRUE.')

    ! trans_coefs = .true., trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept_T, Ypred, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Extended COEFS: TRANS_COEFS=.TRUE., TRANS_X=.TRUE.')

    ! trans_coefs = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept_T, Ypred_T, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Extended COEFS: TRANS_COEFS=.TRUE., TRANS_Y=.TRUE.')

    ! trans_coefs = .true., trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept_T, Ypred_T, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Extended COEFS: TRANS_COEFS=.TRUE., TRANS_X=.TRUE., TRANS_Y=.TRUE.')


end subroutine



subroutine test_pcr_ols_equiv_1d (tests)
    !*  Test models where OLS and PCR yield equivalent results
    class (test_suite) :: tests

    type (test_case), pointer :: tc
    type (lm_config) :: conf_pcr, conf_ols
    type (lm_result) :: res_ols, res_pcr
    real (PREC), dimension(:,:), allocatable :: X
    real (PREC), dimension(:), allocatable :: y
    real (PREC), dimension(:), allocatable :: coefs_true, coefs_ols, coefs_pcr
    real (PREC) ::  intercept_true, intercept_ols, intercept_pcr
    integer :: Nrhs, Nconst, Nobs
    type (status_t) :: status
    type (str) :: msg
    logical :: values_ok
    real (PREC), parameter :: atol = 1.0e-12_PREC, rtol = 1.0e-4_PREC
    integer :: i

    tc => tests%add_test ('PREDICT[1d]: equivalent problems')

    call set_seed (678)

    ! --- Test constant-only model ---

    call random_sample (nrhs=0, nobs=10, X=X, y=y, coefs=coefs_true, &
        intercept=intercept_true, add_intercept=.true., &
        var_error=1.0_PREC, status=status)

    call ols (X=X, Y=Y, res=res_ols, status=status)
    call pcr (X=X, Y=Y, res=res_pcr, status=status)

    values_ok = all_close (res_pcr%intercept, res_ols%intercept, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Constant only, RES present')

    deallocate (X, Y)

    ! --- Test for various number of RHS vars ---

    Nobs = 13

    do i = 1, 10
        Nrhs = i
        call random_sample (nrhs=Nrhs, nobs=Nobs+Nrhs, X=X, y=y, coefs=coefs_true, &
            intercept=intercept_true, add_intercept=.true., &
            var_error=0.1_PREC, status=status)

        ! Force using all available components to get equivalence with OLS
        conf_pcr = lm_config (var_rhs_min=1.0_PREC)

        call ols (X=X, y=y, res=res_ols, status=status)
        status = NF_STATUS_UNDEFINED
        call lm_result_reset (res_pcr)
        call pcr (conf_pcr, X=X, y=y, res=res_pcr, status=status)

        values_ok = all_close (res_pcr%intercept, res_ols%intercept, atol=atol, rtol=rtol) &
            .and. all_close (res_pcr%coefs, res_ols%coefs, atol=atol, rtol=rtol)
        msg = 'Nrhs = ' // str(Nrhs) // '; RES present'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with COEFS, INTERCEPT arguments
        call cond_alloc (coefs_ols, Nrhs)
        call ols (conf_pcr, X=X, y=y, coefs=coefs_ols, intercept=intercept_ols, &
            status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs)
        coefs_pcr(:) = huge(0.0_PREC)
        intercept_pcr = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, intercept=intercept_pcr, &
            status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        msg = 'Nrhs = ' // str(Nrhs) // '; COEFS, INTERCEPT present'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with extended COEFS, INTERCEPT arguments
        Nconst = 1
        call cond_alloc (coefs_ols, Nrhs+Nconst)
        call ols (X=X, y=y, coefs=coefs_ols, intercept=intercept_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs+Nconst)
        coefs_pcr(:) = huge(0.0_PREC)
        intercept_pcr = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, intercept=intercept_pcr, &
            status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        msg = 'Nrhs = ' // str(Nrhs) // '; extended COEFS, INTERCEPT present'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with extended COEFS, no INTERCEPT argument
        Nconst = 1
        call cond_alloc (coefs_ols, Nrhs+Nconst)
        call ols (X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs+Nconst)
        coefs_pcr(:) = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        msg = 'Nrhs = ' // str(Nrhs) // '; extended COEFS, INTERCEPT not present'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with no ADDED intercept
        call cond_alloc (coefs_ols, Nrhs)
        conf_ols = lm_config (add_intercept=.false.)
        call ols (conf_ols, X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs)
        coefs_pcr(:) = huge(0.0_PREC)
        conf_pcr = lm_config (var_rhs_min=1.0_PREC, add_intercept=.false., center=.false.)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (coefs_pcr, coefs_ols, atol=atol, rtol=rtol)
        msg = 'Nrhs = ' // str(Nrhs) // '; no ADDED intercept'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with no added intercept, but intercept in X
        Nconst = 1
        call random_sample (nrhs=Nrhs, nobs=Nobs+Nrhs, X=X, y=y, coefs=coefs_true, &
            add_intercept_x=.true., var_error=0.1_PREC, status=status)

        call cond_alloc (coefs_ols, Nrhs + Nconst)
        conf_ols = lm_config (add_intercept=.false.)
        call ols (conf_ols, X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs + Nconst)
        coefs_pcr(:) = huge(0.0_PREC)
        conf_pcr = lm_config (var_rhs_min=1.0_PREC, add_intercept=.false., center=.false.)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (coefs_pcr, coefs_ols, atol=atol, rtol=rtol)
        msg = 'Nrhs = ' // str(Nrhs) // '; intercept included in X'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

    end do

end subroutine



subroutine test_pcr_ols_equiv_2d (tests)
    !*  Test models where OLS and PCR yield equivalent results
    class (test_suite) :: tests

    type (test_case), pointer :: tc
    type (lm_config) :: conf_pcr, conf_ols
    type (lm_result) :: res_ols, res_pcr
    real (PREC), dimension(:,:), allocatable :: X, Y
    real (PREC), dimension(:,:), allocatable :: coefs_true, coefs_ols, coefs_pcr
    real (PREC), dimension(:), allocatable ::  intercept_true, intercept_ols, intercept_pcr
    integer :: Nobs, Nrhs, Nlhs, Nconst
    type (status_t) :: status
    type (str) :: msg
    logical :: values_ok
    real (PREC), parameter :: atol = 1.0e-12_PREC, rtol = 1.0e-4_PREC
    integer :: i

    tc => tests%add_test ('PREDICT[1d]: equivalent problems')

    call set_seed (678)

    ! --- Test constant-only model ---

    call random_sample (nrhs=0, nobs=10, Nlhs=3, X=X, y=y, coefs=coefs_true, &
        intercept=intercept_true, add_intercept=.true., &
        var_error=1.0_PREC, status=status)

    call ols (X=X, Y=Y, res=res_ols, status=status)
    status = NF_STATUS_OK
    call pcr (X=X, Y=Y, res=res_pcr, status=status)

    values_ok = all_close (res_pcr%intercept, res_ols%intercept, atol=atol, rtol=rtol)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Constant only, RES present')

    deallocate (X, Y)

    ! --- Test for various number of RHS vars ---

    do i = 1, 10, 3
        Nrhs = i
        Nlhs = modulo (i, 3) + 1
        Nobs = 43 + 2*i
        call random_sample (nrhs=Nrhs, nobs=Nobs, Nlhs=Nlhs, X=X, y=y, &
            coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
            var_error=0.1_PREC, status=status)


        call cond_alloc (intercept_ols, Nlhs)
        call cond_alloc (intercept_pcr, Nlhs)

        ! Force using all available components to get equivalence with OLS
        conf_pcr = lm_config (var_rhs_min=1.0_PREC)

        call ols (X=X, y=y, res=res_ols, status=status)
        status = NF_STATUS_UNDEFINED
        call lm_result_reset (res_pcr)
        call pcr (conf_pcr, X=X, y=y, res=res_pcr, status=status)

        values_ok = all_close (res_pcr%intercept_multi, res_ols%intercept_multi, &
                atol=atol, rtol=rtol) &
            .and. all_close (res_pcr%coefs_multi, res_ols%coefs_multi, atol=atol, rtol=rtol)
        msg = 'Nrhs = ' // str(Nrhs) // '; RES present'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with COEFS, INTERCEPT arguments
        call cond_alloc (coefs_ols, [Nrhs,Nlhs])
        call ols (conf_pcr, X=X, y=y, coefs=coefs_ols, intercept=intercept_ols, &
            status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        intercept_pcr(:) = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, intercept=intercept_pcr, &
            status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        msg = 'Nrhs = ' // str(Nrhs) // '; COEFS, INTERCEPT present'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with extended COEFS, INTERCEPT arguments
        Nconst = 1
        call cond_alloc (coefs_ols, [Nrhs+Nconst,Nlhs])
        call ols (X=X, y=y, coefs=coefs_ols, intercept=intercept_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs+Nconst,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        intercept_pcr(:) = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, intercept=intercept_pcr, &
            status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        msg = 'Nrhs = ' // str(Nrhs) // '; extended COEFS, INTERCEPT present'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with extended COEFS, no INTERCEPT argument
        Nconst = 1
        call cond_alloc (coefs_ols, [Nrhs+Nconst,Nlhs])
        call ols (X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs+Nconst,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        msg = 'Nrhs = ' // str(Nrhs) // '; extended COEFS, INTERCEPT not present'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with no ADDED intercept
        call cond_alloc (coefs_ols, [Nrhs,Nlhs])
        conf_ols = lm_config (add_intercept=.false.)
        call ols (conf_ols, X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        conf_pcr = lm_config (var_rhs_min=1.0_PREC, add_intercept=.false., center=.false.)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (coefs_pcr, coefs_ols, atol=atol, rtol=rtol)
        msg = 'Nrhs = ' // str(Nrhs) // '; no ADDED intercept'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        ! Check with no added intercept, but intercept in X
        Nconst = 1
        call random_sample (nrhs=Nrhs, nobs=Nobs+Nrhs, nlhs=Nlhs, X=X, y=y, &
            coefs=coefs_true, add_intercept_x=.true., var_error=0.1_PREC, &
            status=status)

        call cond_alloc (coefs_ols, [Nrhs + Nconst,Nlhs])
        conf_ols = lm_config (add_intercept=.false.)
        call ols (conf_ols, X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs + Nconst,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        conf_pcr = lm_config (var_rhs_min=1.0_PREC, add_intercept=.false., center=.false.)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (coefs_pcr, coefs_ols, atol=atol, rtol=rtol)
        msg = 'Nrhs = ' // str(Nrhs) // '; intercept included in X'
        call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

    end do

end subroutine



subroutine test_pcr_equiv_1d (tests)
    !*  Test variations of equivalent calls to PCR
    class (test_suite) :: tests

    type (test_case), pointer :: tc
    real (PREC), dimension(:,:), allocatable :: X, X_T
    real (PREC), dimension(:), allocatable :: coefs, coefs_ok, coefs_true
    real (PREC), dimension(:), allocatable :: coefs_intercept_ok, coefs_intercept
    real (PREC), dimension(:), allocatable :: y
    real (PREC) :: intercept, intercept_ok, intercept_true
    integer :: Nrhs, Nobs, Ncomp
    real (PREC), parameter :: atol = 1.0e-12_PREC, rtol=1.0e-4_PREC
    logical :: values_ok
    type (status_t) :: status
    type (lm_result) :: res_ok, res
    type (lm_config) :: conf
    type (str) :: msg

    tc => tests%add_test ('PCR[1d] equivalence tests')

    call set_seed (2345)

    do Nrhs = 1, 5

        Nobs = Nrhs + 4
        call random_sample (nrhs=Nrhs, nobs=Nobs, X=X, y=y, &
            coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
            var_error=0.1_PREC, status=status)

        call cond_alloc (X_T, [Nrhs, Nobs])
        X_T(:,:) = transpose (X)

        call cond_alloc (coefs_ok, size(coefs_true))
        call cond_alloc (coefs, size(coefs_true))
        call cond_alloc (coefs_intercept_ok, size(coefs_true) + 1)
        call cond_alloc (coefs_intercept, size(coefs_true) + 1)

        do Ncomp = 0, Nrhs
            conf = lm_config (ncomp=Ncomp)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X, y, coefs=coefs_ok, intercept=intercept_ok, status=status)
            call pcr (conf, X, y, coefs=coefs_intercept_ok, status=status)
            call pcr (conf, X, y, res=res_ok, status=status)

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_x=.true.)
            call pcr (conf, X_T, y, coefs=coefs, intercept=intercept, status=status)

            values_ok = all_close (coefs, coefs_ok, atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; COEFS, INTERCEPT present, TRANS_X=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            status = NF_STATUS_UNDEFINED
            call lm_result_reset (res)
            call pcr (conf, X_T, y, res=res, status=status)

            values_ok = all_close (res%coefs, res_ok%coefs, atol=atol, rtol=rtol) .and. &
                all_close (res%intercept, res_ok%intercept, atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; RES present, TRANS_X=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X_T, y, coefs=coefs_intercept, status=status)

            values_ok = all_close (coefs_intercept, coefs_intercept_ok, &
                atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; extended COEFS, TRANS_X=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)
        end do
    end do

end subroutine



subroutine test_pcr_equiv_2d (tests)
    !*  Test variations of equivalent calls to PCR
    class (test_suite) :: tests

    type (test_case), pointer :: tc
    real (PREC), dimension(:,:), allocatable :: X, X_T
    real (PREC), dimension(:,:), allocatable :: coefs, coefs_T, coefs_ok, coefs_true
    real (PREC), dimension(:,:), allocatable :: coefs_intercept_ok, &
        coefs_intercept, coefs_intercept_T
    real (PREC), dimension(:,:), allocatable :: Y, Y_T
    real (PREC), dimension(:), allocatable :: intercept, intercept_ok, intercept_true
    integer :: Nrhs, Nobs, Ncomp, Nlhs
    real (PREC), parameter :: atol = 1.0e-12_PREC, rtol=1.0e-4_PREC
    logical :: values_ok
    type (status_t) :: status
    type (lm_result) :: res_ok, res
    type (lm_config) :: conf
    type (str) :: msg

    tc => tests%add_test ('PCR[2d] equivalence tests')

    call set_seed (2345)

    do Nrhs = 1, 15, 3

        Nlhs = modulo (Nrhs, 5) + 1

        Nobs = Nrhs + 4
        call random_sample (nrhs=Nrhs, nobs=Nobs, nlhs=Nlhs, X=X, y=y, &
            coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
            var_error=0.1_PREC, status=status)

        call cond_alloc (X_T, [Nrhs, Nobs])
        X_T(:,:) = transpose (X)

        call cond_alloc (Y_T, [Nlhs, Nobs])
        Y_T(:,:) = transpose (Y)

        call cond_alloc (coefs_ok, [Nrhs, Nlhs])
        call cond_alloc (coefs, [Nrhs, Nlhs])
        call cond_alloc (coefs_intercept_ok, [Nrhs + 1, Nlhs])
        call cond_alloc (coefs_intercept, [Nrhs + 1, Nlhs])
        call cond_alloc (coefs_T, [Nlhs, Nrhs])
        call cond_alloc (coefs_intercept_T, [Nlhs, Nrhs + 1])
        call cond_alloc (intercept_ok, Nlhs)
        call cond_alloc (intercept, Nlhs)

        do Ncomp = 0, Nrhs
            conf = lm_config (ncomp=Ncomp)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X, y, coefs=coefs_ok, intercept=intercept_ok, status=status)
            call pcr (conf, X, y, coefs=coefs_intercept_ok, status=status)
            call pcr (conf, X, y, res=res_ok, status=status)

            ! --- Transpose X ---

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_x=.true.)
            call pcr (conf, X_T, y, coefs=coefs, intercept=intercept, status=status)

            values_ok = all_close (coefs, coefs_ok, atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; COEFS, INTERCEPT present, TRANS_X=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            status = NF_STATUS_UNDEFINED
            call lm_result_reset (res)
            call pcr (conf, X_T, y, res=res, status=status)

            values_ok = all_close (res%coefs_multi, res_ok%coefs_multi, &
                    atol=atol, rtol=rtol) .and. &
                all_close (res%intercept_multi, res_ok%intercept_multi, &
                    atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; RES present, TRANS_X=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X_T, y, coefs=coefs_intercept, status=status)

            values_ok = all_close (coefs_intercept, coefs_intercept_ok, &
                atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; extended COEFS, TRANS_X=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            ! --- Transpose Y ---

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_y=.true.)
            call pcr (conf, X, Y_T, coefs=coefs, intercept=intercept, status=status)

            values_ok = all_close (coefs, coefs_ok, atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; COEFS, INTERCEPT present, TRANS_Y=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            status = NF_STATUS_UNDEFINED
            call lm_result_reset (res)
            call pcr (conf, X, Y_T, res=res, status=status)

            values_ok = all_close (res%coefs_multi, res_ok%coefs_multi, &
                atol=atol, rtol=rtol) .and. &
                all_close (res%intercept_multi, res_ok%intercept_multi, &
                    atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; RES present, TRANS_Y=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X, Y_T, coefs=coefs_intercept, status=status)

            values_ok = all_close (coefs_intercept, coefs_intercept_ok, &
                atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; extended COEFS, TRANS_Y=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            ! --- Transpose X and Y ---

            status = NF_STATUS_UNDEFINED
            call lm_result_reset (res)
            conf = lm_config (ncomp=Ncomp, trans_x=.true., trans_y=.true.)
            call pcr (conf, X_T, Y_T, res=res, status=status)

            values_ok = all_close (res%coefs_multi, res_ok%coefs_multi, &
                atol=atol, rtol=rtol) .and. &
                all_close (res%intercept_multi, res_ok%intercept_multi, &
                    atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; RES present, TRANS_X=.TRUE., TRANS_Y=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            ! --- Transpose COEFS ---

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_coefs=.true.)
            call pcr (conf, X, Y, coefs=coefs_T, intercept=intercept, status=status)

            values_ok = all_close (transpose(coefs_T), coefs_ok, &
                    atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; COEFS, INTERCEPT present, TRANS_COEFS=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X, Y, coefs=coefs_intercept_T, status=status)

            values_ok = all_close (transpose(coefs_intercept_T), &
                coefs_intercept_ok, atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; extended COEFS, TRANS_COEFS=.TRUE.'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            ! --- Transpose X, Y and COEFS ---

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_x=.true., trans_y=.true., trans_coefs=.true.)
            call pcr (conf, X_T, Y_T, coefs=coefs_T, intercept=intercept, status=status)

            values_ok = all_close (transpose(coefs_T), coefs_ok, &
                atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; COEFS, INTERCEPT present; X, Y and COEFS transposed'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X_T, Y_T, coefs=coefs_intercept_T, status=status)

            values_ok = all_close (transpose(coefs_intercept_T), &
                coefs_intercept_ok, atol=atol, rtol=rtol)

            msg = 'Nrhs = ' // str(Nrhs) // ', Ncomp = ' // str(Ncomp) &
                // '; extended COEFS; X, Y and COEFS transposed'
            call tc%assert_true (status == NF_STATUS_OK .and. values_ok, msg)

        end do
    end do

end subroutine

end
