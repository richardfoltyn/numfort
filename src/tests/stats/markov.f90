
program test_stats_markov

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_stats
    use numfort_common_testing
    use numfort_common_workspace, workspace => workspace_real64

    use fcore_strings
    use fcore_testing

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("stats_markov unit tests")

    call test_tauchen (tests)
    call test_rouwenhorst_iid (tests)
    call test_rouwenhorst_python (tests)

    call tests%print ()

end subroutine


subroutine test_tauchen (tests)
    !*  TEST_TAUCHEN compares discretized moments of Markov chain generated
    !   using the TAUCHEN method to those tabulated in Tauchen (1986).
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC) :: rho, sigma
    real (PREC) :: rho_m, sigma_m
    real (PREC), dimension(:), allocatable :: states
    real (PREC), dimension(:,:), allocatable :: tm
    integer :: n
    type (status_t) :: status
    
    tc => tests%add_test ('TAUCHEN unit tests')

    ! Replicate results reported in Tauchen (1986) paper
    rho = 0.10d0
    sigma = 0.101
    n = 9

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, sigma_cond=.false., status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    call tc%assert_true (all_close (rho_m, 0.1d0, atol=1.0d-3) .and. &
        all_close (sigma_m, 0.103d0, atol=1.0d-3), &
        "Moments for rho=" // str(rho, 'f4.2') // ", sigma=" // str(sigma, 'f5.3'))
    deallocate (states, tm)

    rho = 0.80d0
    sigma = 0.167
    n = 9

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, sigma_cond=.false., status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    call tc%assert_true (all_close (rho_m, 0.798d0, atol=1.0d-3) .and. &
        all_close (sigma_m, 0.176d0, atol=1.0d-3), &
        "Moments for rho=" // str(rho, 'f4.2') // ", sigma=" // str(sigma, 'f5.3'))
    deallocate (states, tm)

    rho = 0.90d0
    sigma = 0.229
    n = 9

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, sigma_cond=.false., status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    call tc%assert_true (all_close (rho_m, 0.898d0, atol=1.0d-3) .and. &
        all_close (sigma_m, 0.253d0, atol=1.0d-3), &
        "Moments for rho=" // str(rho, 'f4.2') // ", sigma=" // str(sigma, 'f5.3'))
    deallocate (states, tm)

    rho = 0.90d0
    sigma = 0.229
    n = 5

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, sigma_cond=.false., status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    call tc%assert_true (all_close (rho_m, 0.932d0, atol=1.0d-3) .and. &
        all_close (sigma_m, 0.291d0, atol=1.0d-3), &
        "Moments for rho=" // str(rho, 'f4.2') // ", sigma=" // str(sigma, 'f5.3'))
    deallocate (states, tm)

    ! Test some degenerate processes
    rho = 0.0d0
    sigma = sqrt(0.0522d0)
    n = 5

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    call tc%assert_true (all_close (rho_m, 0.0d0, atol=1.0d-8), &
        "Markov chain with rho=0")
    deallocate (states, tm)

end subroutine


subroutine test_rouwenhorst_iid (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC) :: rho, sigma
    real (PREC) :: rho_m, sigma_m
    real (PREC), dimension(:), allocatable :: states
    real (PREC), dimension(:,:), allocatable :: tm
    integer :: n
    type (status_t) :: status

    tc => tests%add_test ('ROUWENHORST (IID) unit tests')

    rho = 0.0d0
    sigma = sqrt(0.0522d0)
    n = 5
    allocate (states(n), tm(n,n))
    call rouwenhorst (rho, sigma, states, tm, status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    call tc%assert_true (all_close (rho_m, 0.0d0, atol=1.0d-8), &
        "Markov chain with rho=0")
    deallocate (states, tm)


end subroutine


subroutine test_rouwenhorst_python (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC) :: rho, sigma
    real (PREC), dimension(:), allocatable :: states, edist, edist2
    real (PREC), dimension(:,:), allocatable :: tm, tm_T
    type (workspace) :: work
    integer :: n
    type (status_t) :: status

    integer, parameter :: N2 = 9
    real (PREC), parameter :: TM2_desired(*,*) = reshape([ &
        6.63420431d-01,   2.79334918d-01,   5.14564323d-02,   5.41646656d-03, &
        3.56346484d-04,   1.50040625d-05,   3.94843750d-07,   5.93750000d-09, &
        3.90625000d-11,   3.49168648d-02,   6.76284539d-01,   2.46449229d-01, &
        3.87704975d-02,   3.39466914d-03,   1.78469375d-04,   5.63171875d-06, &
        9.87500000d-08,   7.42187500d-10,   1.83772973d-03,   7.04140653d-02, &
        6.85549548d-01,   2.12408226d-01,   2.77697840d-02,   1.94249469d-03, &
        7.65292188d-05,   1.60906250d-06,   1.41015625d-08,   9.67226172d-05, &
        5.53864250d-03,   1.06204113d-01,   6.91139238d-01,   1.77494044d-01, &
        1.85302288d-02,   9.71247344d-04,   2.54956250d-05,   2.67929688d-07, &
        5.09066406d-06,   3.87962188d-04,   1.11079136d-02,   1.41995235d-01, &
        6.93007596d-01,   1.41995235d-01,   1.11079136d-02,   3.87962188d-04, &
        5.09066406d-06,   2.67929688d-07,   2.54956250d-05,   9.71247344d-04, &
        1.85302288d-02,   1.77494044d-01,   6.91139238d-01,   1.06204113d-01, &
        5.53864250d-03,   9.67226172d-05,   1.41015625d-08,   1.60906250d-06, &
        7.65292188d-05,   1.94249469d-03,   2.77697840d-02,   2.12408226d-01, &
        6.85549548d-01,   7.04140653d-02,   1.83772973d-03,   7.42187500d-10, &
        9.87500000d-08,   5.63171875d-06,   1.78469375d-04,   3.39466914d-03, &
        3.87704975d-02,   2.46449229d-01,   6.76284539d-01,   3.49168648d-02, &
        3.90625000d-11,   5.93750000d-09,   3.94843750d-07,   1.50040625d-05, &
        3.56346484d-04,   5.41646656d-03,   5.14564323d-02,   2.79334918d-01, &
        6.63420431d-01], &
        shape=[N2,N2], order=[2, 1])

    real (PREC), parameter :: states2_desired(*) = [ &
        -0.64888568d0, -0.48666426d0, -0.32444284d0, -0.16222142d0, 0.0d0,  &
         0.16222142d0,  0.32444284d0,  0.48666426d0,  0.64888568d0]

    real (PREC), parameter :: EDIST2_DESIRED(*) = [ &
        0.00390625d0, 0.03125d0, 0.109375d0, 0.21875d0, 0.2734375d0, &
        0.21875d0, 0.109375d0, 0.03125d0, 0.00390625d0 ]

    tc => tests%add_test ('ROUWENHORST Fortran vs. Python unit tests')

    ! Compare to results from Python implementation for N=9, rho=0.9, sigma=0.1
    rho = 0.9d0
    sigma = 0.1d0
    n = N2

    allocate (states(n), tm(n,n), tm_T(n,n), edist(n), edist2(n))

    call rouwenhorst (rho, sigma, states, tm, status=status)

    call tc%assert_true (status == NF_STATUS_OK, &
        "ROUWENHORST exit code")

    call tc%assert_true (all_close (tm, TM2_desired, rtol=1.0d-7), &
        "TM close to Python result")

    call tc%assert_true (all_close (states, states2_desired, atol=1.0d-7), &
        "States close to Python result")

    call markov_ergodic_dist (tm, edist, inverse=.true., work=work, status=status)

    call tc%assert_true (status == NF_STATUS_OK, &
        "ERGODIC_DIST exit code")

    call tc%assert_true (all_close (edist, EDIST2_DESIRED, atol=1.0d-4), &
        "Ergodic dist. close to Python result")

    tm_T(:,:) = transpose(tm)
    status = NF_STATUS_UNDEFINED
    call markov_ergodic_dist (tm_T, edist2, inverse=.true., is_transposed=.true., &
        work=work, status=status)

    call tc%assert_true (status == NF_STATUS_OK, &
        "ERGODIC_DIST (is_transposed=.true.) exit code")

    call tc%assert_true (all_close (edist2, edist, atol=1.0d-10, rtol=0.0d0), &
        "ERGODIC_DIST (is_transposed=.true.) result")

    edist2(:) = 0.0
    status = NF_STATUS_UNDEFINED
    call markov_ergodic_dist (tm, edist2, inverse=.false., work=work, status=status)

    call tc%assert_true (status == NF_STATUS_OK, &
        "ERGODIC_DIST (inverse=.false.) exit code")

    call tc%assert_true (all_close (edist, edist, atol=1.0d-10, rtol=0.0d0), &
        "ERGODIC_DIST (inverse=.false.) result")

    edist2(:) = 0.0
    status = NF_STATUS_UNDEFINED
    call markov_ergodic_dist (tm_T, edist2, inverse=.false., is_transposed=.true., &
        work=work, status=status)

    call tc%assert_true (status == NF_STATUS_OK, &
        "ERGODIC_DIST (inverse=.false., is_transposed=.true.) exit code")

    call tc%assert_true (all_close (edist, edist, atol=1.0d-10, rtol=0.0d0), &
        "ERGODIC_DIST (inverse=.false., is_transposed=.true.) result")

    deallocate (states, tm, tm_T, edist, edist2)


end subroutine

end program
