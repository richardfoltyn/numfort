
program test_stats_markov

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_stats
    use numfort_common_testing

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
    call test_rouwenhorst (tests)

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


subroutine test_rouwenhorst (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC) :: rho, sigma
    real (PREC) :: rho_m, sigma_m
    real (PREC), dimension(:), allocatable :: states
    real (PREC), dimension(:,:), allocatable :: tm
    integer :: n
    type (status_t) :: status

    tc => tests%add_test ('ROUWENHORST unit tests')

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

end program
