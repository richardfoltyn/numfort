

program test_stats_lm
    !*  Unit tests for OLS estimation routines in numfort_stats_lm module.

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_testing
    use numfort_stats, only: set_seed
    use numfort_stats_lm, lm_data => lm_data_real64
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

    call test_ridge_1d_args (tests)
    call test_ridge_2d_args (tests)
    call test_ridge_1d (tests)

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
    type (status_t) :: status

    tc => tests%add_test ("OLS[1d] argument checking")

    ! === incompatible X, Y arguments ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs+1), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, status=status)
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
    call ols (y, x, beta, trans_x=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, Y")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs-1))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs+1,nobs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, trans_x=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS with ADD_CONST=.FALSE. ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, add_const=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS with ADD_CONST=.FALSE.")

    deallocate (x, y, beta)

    ! === NOBS < NCOEFS ===
    nrhs = 3
    nobs = nrhs
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, add_const=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "NOBS < NCOEFS")

    deallocate (x, y, beta)

    ! === Invalid RCOND ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, rcond=-0.1_PREC, status=status)
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
    call ols (y, x, beta, status=status)
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
    call ols (y, x, beta, trans_x=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, Y")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs-1,nlhs))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs+1,nobs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, trans_x=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS with ADD_CONST=.FALSE. ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, add_const=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS with ADD_CONST=.FALSE.")

    deallocate (x, y, beta)

    ! === incompabible Y, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs-1))

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, status=status)
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
    call ols (y, x, beta, status=status)
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
    call ols (y, x, beta, rcond=-0.1_PREC, status=status)
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
    real (PREC), parameter :: RTOL = 1.0e-3_PREC, ATOL=0.0_PREC
        ! Arguments for ALL_CLOSE()
    real (PREC) :: rsq
    logical :: all_ok
    type (lm_data), dimension(:), allocatable :: lm
    type (status_t) :: status

    tc => tests%add_test ('OLS[2d] unit tests')

    call set_seed (123)

    ! === exactly determined system with NOBS = NRHS, no const ===
    nrhs = 10
    nobs = 10
    ncoefs = nrhs
    nlhs = 3

    allocate (x(nobs,nrhs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))
    allocate (lm(nlhs))

    call random_number (x)
    call random_number (beta_ok)

    call BLAS95_GEMM (x, beta_ok, y)

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, add_const=.false., rank=rank, res=lm, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        'Exactly identified system: estimation')

    ! Test without optional COEFS argument
    deallocate (lm)
    allocate (lm(nlhs))
    status = NF_STATUS_UNDEFINED
    call ols (y, x, add_const=.false., rank=rank, res=lm, status=status)
    all_ok = .true.
    do i = 1, nlhs
        all_ok = all_ok .and. all_close (lm(i)%coefs, beta_ok(:,i), rtol=RTOL, atol=ATOL)
    end do
    call tc%assert_true (status == NF_STATUS_OK .and. all_ok, &
        'Exactly identified system: estimation w/o COEFS argument')

    ! Compute R^2 which must be 1.0 in this case
    all_ok = .true.
    do i = 1, nlhs
        status = NF_STATUS_UNDEFINED
        call post_estim (lm(i), y(:,i), x, rsq, status)
        all_ok = all_ok .and. abs(rsq-1.0_PREC) < 1.0e-12_PREC &
            .and. status == NF_STATUS_OK
    end do
    call tc%assert_true (all_ok, 'Exactly identified system: R^2')

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
    call rvs (norm, y, loc=0.0_PREC, scale=1.0e-3_PREC)

    ! Create LHS variables
    call BLAS95_GEMM (x, beta_ok, y, beta=1.0_PREC)

    status = NF_STATUS_OK
    call ols (y, x, beta, add_const=.false., rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        "Underdetermined system with small residuals")

    deallocate (x, y, beta, beta_ok)


    ! === overdetermined system with intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs + 1
    nlhs = 11

    allocate (x(nobs,nrhs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)
    call rvs (norm, y, loc=0.0_PREC, scale=1.0e-3_PREC)

    ! Create LHS variables
    forall (i=1:nlhs) y(:,i) = y(:,i) + beta_ok(1,i)
    allocate (work(nrhs,nlhs), source=beta_ok(2:,:))
    call BLAS95_GEMM (x, work, y, beta=1.0_PREC)
    deallocate (work)

    status = NF_STATUS_OK
    call ols (y, x, beta, add_const=.true., rank=rank, status=status)
    ! Apply less strict check on beta, value of 1.0e-3 does not seem to work
    ! for NLHS=11
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=1.0e-2_PREC, atol=ATOL), &
        "Underdetermined system with small residuals, ADD_CONST=.TRUE.")

    deallocate (x, y, beta, beta_ok)


    ! === overdetermined system with intercept, transposed X ===
    nrhs = 3
    nobs = 500
    ncoefs = nrhs + 1
    nlhs = 7

    allocate (x(nrhs,nobs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)
    call rvs (norm, y, loc=0.0_PREC, scale=1.0e-3_PREC)

    ! Create LHS variables
    forall (i=1:nlhs) y(:,i) = y(:,i) + beta_ok(1,i)
    allocate (work(nrhs,nlhs), source=beta_ok(2:,:))
    call BLAS95_GEMM (x, work, y, beta=1.0_PREC, transa='T')
    deallocate (work)

    status = NF_STATUS_OK
    call ols (y, x, beta, add_const=.true., trans_x=.true., rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=1.0e-2_PREC, atol=ATOL), &
        "Underdetermined system with small residuals, ADD_CONST=.TRUE., TRANS_X=.TRUE.")

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
    real (PREC), parameter :: RTOL = 1.0e-3_PREC, ATOL=0.0_PREC
        ! Arguments for ALL_CLOSE()
    real (PREC) :: rsq
    type (status_t) :: status
    type (lm_data), allocatable :: lm

    tc => tests%add_test ('OLS[1d] unit tests')

    call set_seed (456)

    ! === exactly determined system with NOBS = NRHS, no const ===
    nrhs = 10
    nobs = 10
    ncoefs = nrhs

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_ok(ncoefs))

    call random_number (x)
    call random_number (beta_ok)

    call BLAS95_GEMV (x, beta_ok, y)

    status = NF_STATUS_UNDEFINED
    allocate (lm)
    call ols (y, x, beta, add_const=.false., rank=rank, res=lm, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        'Exactly identified system: estimation')

    ! Test w/o optional COEFS argument, coefs are stored in LM_DATA
    deallocate (lm)
    allocate (lm)
    status = NF_STATUS_UNDEFINED
    call ols (y, x, add_const=.false., rank=rank, res=lm, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (lm%coefs, beta_ok, rtol=RTOL, atol=ATOL), &
        'Exactly identifies system: estimation w/o COEFS argument')

    ! Compute R^2 which must be 1.0 in this case
    status = NF_STATUS_UNDEFINED
    call post_estim (lm, y, x, rsq, status)
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
    call rvs (norm, y, loc=0.0_PREC, scale=1.0e-3_PREC)

    ! Create LHS variables
    call BLAS95_GEMV (x, beta_ok, y, beta=1.0_PREC)

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, add_const=.false., rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        "Underdetermined system with small residuals")

    deallocate (x, y, beta, beta_ok)


    ! === overdetermined system with intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs + 1

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_ok(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)
    call rvs (norm, y, loc=0.0_PREC, scale=1.0e-3_PREC)

    ! Create LHS variables
    y(:) = y(:) + beta_ok(1)
    allocate (work(nrhs), source=beta_ok(2:))
    call BLAS95_GEMV (x, work, y, beta=1.0_PREC)
    deallocate (work)

    status = NF_STATUS_UNDEFINED
    call ols (y, x, beta, add_const=.true., rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        "Underdetermined system with small residuals, ADD_CONST=.TRUE.")

    deallocate (x, y, beta, beta_ok)


    ! === overdetermined system with intercept, transposed X ===
    nrhs = 2
    nobs = 500
    ncoefs = nrhs + 1

    allocate (x(nrhs,nobs), y(nobs), beta(ncoefs), beta_ok(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)
    call rvs (norm, y, loc=0.0_PREC, scale=1.0e-3_PREC)

    ! Create LHS variables
    y(:) = y(:) + beta_ok(1)
    allocate (work(nrhs), source=beta_ok(2:))
    call BLAS95_GEMV (x, work, y, beta=1.0_PREC, trans='T')
    deallocate (work)

    status = NF_STATUS_OK
    call ols (y, x, beta, add_const=.true., trans_x=.true., rank=rank, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (beta, beta_ok, rtol=RTOL, atol=ATOL), &
        "Underdetermined system with small residuals, ADD_CONST=.TRUE., TRANS_X=.TRUE.")

    deallocate (x, y, beta, beta_ok)

end subroutine



subroutine test_ridge_1d_args (tests)
    !*  Unit tests for rigde regression argument checking
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: x
    real (PREC), dimension(:), allocatable :: y
    real (PREC), dimension(:), allocatable :: beta

    real (PREC) :: lambda
    integer :: nrhs, nobs, ncoefs
    type (status_t) :: status

    tc => tests%add_test ("Ridge[1d] argument checking")

    lambda = 0.0

    ! === incompatible X, Y arguments ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs+1), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, status=status)
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
    call ridge (y, x, lambda, beta, trans_x=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, Y")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs-1))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs+1,nobs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, trans_x=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS with ADD_CONST=.FALSE. ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, add_const=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS with ADD_CONST=.FALSE.")

    deallocate (x, y, beta)

    ! === NOBS < NCOEFS ===
    nrhs = 3
    nobs = nrhs
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, add_const=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "NOBS < NCOEFS")

    deallocate (x, y, beta)

    ! === Invalid RCOND ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, rcond=-0.1_PREC, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Invalid RCOND")

    deallocate (x, y, beta)

    ! === Invalid lambda ===

    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    lambda = -1.0_PREC
    call ridge (y, x, lambda, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Invalid LAMBDA")

    deallocate (x, y, beta)


end subroutine


subroutine test_ridge_2d_args (tests)
    !*  Unit tests for ridge regression argument checking
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: x, y
    real (PREC), dimension(:,:), allocatable :: beta

    integer :: nrhs, nlhs, nobs, ncoefs
    real (PREC) :: lambda
    type (status_t) :: status

    tc => tests%add_test ("Ridge[2d] argument checking")

    lambda = 0.0

    ! === incompatible X, Y arguments ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs+1,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, status=status)
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
    call ridge (y, x, lambda, beta, trans_x=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, Y")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs-1,nlhs))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs+1,nobs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, trans_x=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X^T, COEFS")

    deallocate (x, y, beta)

    ! === incompabible X, COEFS with ADD_CONST=.FALSE. ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, add_const=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable arrays X, COEFS with ADD_CONST=.FALSE.")

    deallocate (x, y, beta)

    ! === incompabible Y, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs-1))

    status = NF_STATUS_UNDEFINED
    call ridge (y, x, lambda, beta, status=status)
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
    call ridge (y, x, lambda, beta, status=status)
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
    call ridge (y, x, lambda, beta, rcond=-0.1_PREC, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Invalid RCOND")

    deallocate (x, y, beta)

    ! === Invalid LAMBDA ===
    nrhs = 3
    nlhs = 2
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    lambda = -1.0
    call ridge (y, x, lambda, beta, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Invalid LAMBDA")

    deallocate (x, y, beta)

end subroutine



subroutine test_ridge_1d (tests)
    !*  Unit tests for various OLS regression problems using the 1d API
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: x
    real (PREC), dimension(:), allocatable :: y, beta, beta_ok, work
    real (PREC), dimension(:), allocatable :: beta_ols
    integer :: nrhs, ncoefs, nobs
    integer :: rank, rank_ols
    real (PREC), parameter :: RTOL = 1.0e-3_PREC, ATOL=0.0_PREC
    logical :: add_const, values_ok, trans_x
    ! Arguments for ALL_CLOSE()
    real (PREC) :: rsq, lambda, ssq, ssq_ols
    type (status_t) :: status
    type (lm_data), allocatable :: lm, lm_ols

    tc => tests%add_test ('Ridge[1d] unit tests')

    call set_seed (456)

    ! === exactly determined system with NOBS = NRHS, no const ===

    nrhs = 10
    nobs = 10
    ncoefs = nrhs

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_ok(ncoefs))
    allocate (beta_ols(ncoefs))

    call random_number (x)
    call random_number (beta_ok)

    call BLAS95_GEMV (x, beta_ok, y)

    allocate (lm, lm_ols)
    add_const = .false.

    call ols (y, x, beta_ols, add_const=add_const, rank=rank_ols, res=lm_ols)

    ! Test with lambda = 0.0
    status = NF_STATUS_UNDEFINED
    lambda = 0.0
    call ridge (y, x, lambda, beta, add_const=add_const, rank=rank, res=lm, status=status)

    values_ok = all_close (beta, beta_ok, rtol=RTOL, atol=ATOL)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Estim: Exactly identified system, LAMBDA=0')

    ! Test with lambda > 0.0
    status = NF_STATUS_UNDEFINED
    lambda = 0.1
    beta(:) = 0.0
    call ridge (y, x, lambda, beta, add_const=add_const, rank=rank, res=lm, status=status)

    ssq = sum(beta**2.0)
    ssq_ols = sum(beta_ols**2.0)
    call tc%assert_true (status == NF_STATUS_OK .and. ssq <= ssq_ols, &
        'Estim: Exactly identified system, LAMBDA>0')

    deallocate (x, y, beta, beta_ok, beta_ols)

    ! === overdetermined system without intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_ok(ncoefs))
    allocate (beta_ols(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)
    call rvs (norm, y, loc=0.0_PREC, scale=1.0e-3_PREC)

    ! Create LHS variables
    call BLAS95_GEMV (x, beta_ok, y, beta=1.0_PREC)

    add_const = .false.

    call ols (y, x, beta_ols, add_const=add_const, rank=rank_ols)

    ! Lambda = 0
    status = NF_STATUS_UNDEFINED
    lambda = 0.0
    call ridge (y, x, lambda, beta, add_const=add_const, rank=rank, status=status)

    values_ok = all_close (beta, beta_ok, rtol=RTOL, atol=ATOL)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        "Estim: Underdetermined system with small residuals, LAMBDA=0")

    ! Lambda > 0
    status = NF_STATUS_UNDEFINED
    lambda = 0.1d0
    call ridge (y, x, lambda, beta, add_const=add_const, rank=rank, status=status)

    ssq = sum(beta**2.0)
    ssq_ols = sum(beta_ols**2.0)
    call tc%assert_true (status == NF_STATUS_OK .and. ssq <= ssq_ols, &
        "Estim: Underdetermined system with small residuals, LAMBDA>0")

    deallocate (x, y, beta, beta_ok, beta_ols)


    ! === overdetermined system with intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs + 1

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_ok(ncoefs))
    allocate (beta_ols(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)
    call rvs (norm, y, loc=0.0_PREC, scale=1.0e-3_PREC)

    ! Create LHS variables
    y(:) = y(:) + beta_ok(1)
    allocate (work(nrhs), source=beta_ok(2:))
    call BLAS95_GEMV (x, work, y, beta=1.0_PREC)
    deallocate (work)

    add_const = .true.

    call ols (y, x, beta_ols, add_const=add_const, rank=rank_ols)

    ! Lambda = 0.0
    status = NF_STATUS_UNDEFINED
    lambda = 0.0
    call ridge (y, x, lambda, beta, add_const=add_const, rank=rank, status=status)

    values_ok = all_close (beta, beta_ok, rtol=RTOL, atol=ATOL)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        "Estim: Underdetermined system with small residuals, LAMBDA=0, ADD_CONST=.TRUE.")

    ! Lambda > 0.0
    status = NF_STATUS_UNDEFINED
    lambda = 0.5
    call ridge (y, x, lambda, beta, add_const=add_const, rank=rank, status=status)

    ssq = sum(beta**2.0)
    ssq_ols = sum(beta_ols**2.0)
    call tc%assert_true (status == NF_STATUS_OK .and. ssq <= ssq_ols, &
        "Estim: Underdetermined system with small residuals, LAMBDA>0, ADD_CONST=.TRUE.")

    deallocate (x, y, beta, beta_ok, beta_ols)


    ! === overdetermined system with intercept, transposed X ===
    nrhs = 2
    nobs = 500
    ncoefs = nrhs + 1

    allocate (x(nrhs,nobs), y(nobs), beta(ncoefs), beta_ok(ncoefs))
    allocate (beta_ols(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)
    call rvs (norm, y, loc=0.0_PREC, scale=1.0e-3_PREC)

    ! Create LHS variables
    y(:) = y(:) + beta_ok(1)
    allocate (work(nrhs), source=beta_ok(2:))
    call BLAS95_GEMV (x, work, y, beta=1.0_PREC, trans='T')
    deallocate (work)

    add_const = .true.
    trans_x = .true.

    call ols (y, x, beta_ols, add_const=add_const, trans_x=trans_x, rank=rank_ols)

    ! Lambda = 0.0
    status = NF_STATUS_UNDEFINED
    lambda = 0.0
    call ridge (y, x, lambda, beta, add_const=add_const, trans_x=trans_x, rank=rank, status=status)

    values_ok = all_close (beta, beta_ok, rtol=RTOL, atol=ATOL)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        "Estim: Underdet. system with small resid, LAMBDA=0, ADD_CONST=.TRUE., TRANS_X=.TRUE.")

    ! Lambda > 0.0
    status = NF_STATUS_UNDEFINED
    lambda = 0.1
    call ridge (y, x, lambda, beta, add_const=add_const, trans_x=trans_x, rank=rank, status=status)

    ssq = sum(beta**2.0)
    ssq_ols = sum(beta_ols**2.0)
    call tc%assert_true (status == NF_STATUS_OK .and. ssq <= ssq_ols, &
        "Estim: Underdet. system with small resid, LAMBDA>0, ADD_CONST=.TRUE., TRANS_X=.TRUE.")


    deallocate (x, y, beta, beta_ok, beta_ols)

end subroutine


end
