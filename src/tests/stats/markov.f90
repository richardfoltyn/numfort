
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

    call test_simulate (tests)
    call test_simulate_advanced (tests)

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
    call markov_ergodic_dist (tm_T, edist2, inverse=.true., transposed=.true., &
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
    call markov_ergodic_dist (tm_T, edist2, inverse=.false., transposed=.true., &
        work=work, status=status)

    call tc%assert_true (status == NF_STATUS_OK, &
        "ERGODIC_DIST (inverse=.false., is_transposed=.true.) exit code")

    call tc%assert_true (all_close (edist, edist, atol=1.0d-10, rtol=0.0d0), &
        "ERGODIC_DIST (inverse=.false., is_transposed=.true.) result")

    deallocate (states, tm, tm_T, edist, edist2)

end subroutine



subroutine test_simulate (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: transm, transm_T
    real (PREC), dimension(:), allocatable :: states
    integer, dimension(:), allocatable :: sim, sim_T
    integer :: n, Nobs
    integer :: in, ir, is
    real (PREC) :: rho, sigma
    real (PREC), parameter :: RHO_ALL(*) = [ real (PREC) :: 0.0, 0.2, 0.5, 0.9]
    real (PREC), parameter :: SIGMA_ALL(*) = [ real (PREC) :: 0.01, 0.1, 0.2]
    integer, parameter :: N_ALL(*) = [3, 5, 7]
    logical :: transposed, values_ok
    type (str) :: msg
    type (status_t) :: status, status1, status2

    tc => tests%add_test ('Unit tests for Markov simulation routine')

    ! === Input checks ===

    ! 1. Invalid transition matrix
    Nobs = 100
    allocate (sim(Nobs))
    n = 2
    allocate (transm(n,n), source=0.0_PREC)
    allocate (transm_T(n,n), source=0.0_PREC)

    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 1, sim, transposed=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid transm. #1, transposed=.false.')

    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 1, sim, transposed=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid transm. #1, transposed=.true.')

    transm(:,:) = 1.0
    transm_T(:,:) = 1.0
    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 1, sim, transposed=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid transm. #2, transposed=.false.')

    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 1, sim, transposed=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid transm. #2, transposed=.true.')

    ! 2. Invalid initial value
    transm(:,:) = 1.0_PREC / n
    transm_T(:,:) = 1.0_PREC / n
    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 0, sim, transposed=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid initial state < 0, transposed=.false.')

    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, n+1, sim, transposed=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid initial state > N, transposed=.false.')

    deallocate (transm, transm_T, sim)

    ! === Compare transposed vs. non-transposed TM input ===

    Nobs = 10000
    allocate (sim(Nobs), sim_T(Nobs))

    do in = 1, size(N_ALL)

        n = N_ALL(in)

        allocate (transm(n,n), transm_T(n,n), states(n))

        do ir = 1, size(RHO_ALL)
            do is = 1, size(SIGMA_ALL)
                rho = RHO_ALL(ir)
                sigma = SIGMA_ALL(is)

                call rouwenhorst (rho, sigma, states, transm)
                transm_T(:,:) = transpose(transm)

                call set_seed (1234)
                status1 = NF_STATUS_UNDEFINED
                sim(:) = 0
                call markov_simulate (transm, 1, sim, transposed=.false., &
                    status=status1)

                call set_seed (1234)
                status2 = NF_STATUS_UNDEFINED
                sim_T(:) = 0
                call markov_simulate (transm_T, 1, sim_T, transposed=.true., &
                    status=status2)

                values_ok = all(sim == sim_T)

                msg = 'Simulation for N=' // str(n, 'i0') &
                    // '; rho=' // str(rho, 'f0.2') &
                    // '; sigma=' // str(sigma, 'f0.2') &
                    // ': transposed == non-transposed'
                call tc%assert_true (values_ok .and. status1 == NF_STATUS_OK &
                    .and. status2 == NF_STATUS_OK, msg)

            end do
        end do

        deallocate (transm, transm_T, states)
    end do

end subroutine



subroutine test_simulate_advanced (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: transm, transm_T
    real (PREC), dimension(:,:), allocatable :: pmf_trans, pmf_sample
    real (PREC), dimension(:), allocatable :: states, edist
    integer, dimension(:,:), allocatable :: ntrans
    integer, dimension(:), allocatable :: sim, sim_T
    integer :: n, Nobs
    integer :: in, ir, is, i
    real (PREC) :: rho, sigma
    real (PREC), parameter :: RHO_ALL(*) = [ real (PREC) :: 0.0, 0.2d0, 0.5d0, 0.9d0]
    real (PREC), parameter :: SIGMA_ALL(*) = [ real (PREC) :: 0.01d0, 0.1d0, 0.2d0]
    integer, parameter :: N_ALL(*) = [3, 5]
    real (PREC), parameter :: rtol = 5.0e-2_PREC, atol=1.0e-3_PREC
    logical :: transposed, values_ok
    type (str) :: msg
    type (status_t) :: status, status1, status2

    tc => tests%add_test ('Unit tests for Markov adv. simulation routine')

    ! === Input checks ===

    ! 1. Invalid transition matrix
    Nobs = 100
    allocate (sim(Nobs))
    n = 2
    allocate (transm(n,n), source=0.0_PREC)
    allocate (transm_T(n,n), source=0.0_PREC)

    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 1, sim, transposed=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid transm. #1, transposed=.false.')

    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 1, sim, transposed=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid transm. #1, transposed=.true.')

    transm(:,:) = 1.0
    transm_T(:,:) = 1.0
    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 1, sim, transposed=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid transm. #2, transposed=.false.')

    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 1, sim, transposed=.true., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid transm. #2, transposed=.true.')

    ! 2. Invalid initial value
    transm(:,:) = 1.0_PREC / n
    transm_T(:,:) = 1.0_PREC / n
    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 0, sim, transposed=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid initial state < 0, transposed=.false.')

    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, n+1, sim, transposed=.false., status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        'Invalid initial state > N, transposed=.false.')

    deallocate (transm, transm_T, sim)

    ! === Compare transposed vs. non-transposed TM input ===

    Nobs = 10000
    allocate (sim(Nobs), sim_T(Nobs))

    do in = 1, size(N_ALL)

        n = N_ALL(in)

        allocate (transm(n,n), transm_T(n,n), states(n), edist(n))
        allocate (pmf_trans(n,n), pmf_sample(n,n))
        allocate (ntrans(n,n))

        do ir = 1, size(RHO_ALL)
            do is = 1, size(SIGMA_ALL)
                rho = RHO_ALL(ir)
                sigma = SIGMA_ALL(is)

                call rouwenhorst (rho, sigma, states, transm)
                transm_T(:,:) = transpose(transm)

                call set_seed (1234)
                status1 = NF_STATUS_UNDEFINED
                sim(:) = 0
                call markov_simulate_advanced (transm, 1, sim, rtol=rtol, &
                    atol=atol, transposed=.false., status=status1)

                call set_seed (1234)
                status2 = NF_STATUS_UNDEFINED
                sim_T(:) = 0
                call markov_simulate_advanced (transm_T, 1, sim_T, rtol=rtol, &
                    atol=atol, transposed=.true., status=status2)

                values_ok = all(sim == sim_T)

                msg = 'Simulation for N=' // str(n, 'i0') &
                    // '; rho=' // str(rho, 'f0.2') &
                    // '; sigma=' // str(sigma, 'f0.2') &
                    // ': transposed == non-transposed'
                call tc%assert_true (values_ok .and. status1 == NF_STATUS_OK &
                    .and. status2 == NF_STATUS_OK, msg)

                ! Check that generated sequence is close to actual
                ! transition matrix
                call count_trans (sim, ntrans)

                ! Theoretical transition PMF
                call markov_ergodic_dist (transm, edist, inverse=.true.)
                do i = 1, n
                    pmf_trans(i,:) = transm(i,:) * edist(i)
                end do

                ! Transition PMF for simulated sample
                do i = 1, n
                    pmf_sample(i,:) = ntrans(i,:) / real(sum(ntrans), PREC)
                end do

                ! This is not exactly the terminal condition that simulation
                ! routine computes, so we need to less restrictive tolerance
                values_ok = all_close (pmf_sample, pmf_trans, &
                    atol=atol + 1.0e-3_PREC, rtol=rtol)

                msg = 'Simulation for N=' // str(n, 'i0') &
                        // '; rho=' // str(rho, 'f0.2') &
                        // '; sigma=' // str(sigma, 'f0.2') &
                        // ': close trans. frequencies'
                call tc%assert_true (values_ok, msg)

            end do
        end do

        deallocate (transm, transm_T, ntrans, states, edist)
        deallocate (pmf_trans, pmf_sample)
    end do

end subroutine



pure subroutine count_trans (x, num)
    !*  COUNT_TRANS returns the number of observed transitions for each
    !   transition type.
    integer, intent(in), dimension(:), contiguous :: x
    integer, intent(out), dimension(:,:), contiguous :: num
    !*  Matrix containing the number of transitions for each possible
    !   combination of origin and destrination values.

    integer :: i, n
    integer :: ifrom,  ito

    n = size(x)
    num = 0

    do i = 1, n - 1
        ifrom = x(i)
        ito = x(i+1)
        num(ifrom, ito) = num(ifrom,ito) + 1
    end do

end subroutine



end program
