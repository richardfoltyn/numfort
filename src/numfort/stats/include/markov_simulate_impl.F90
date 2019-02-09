

subroutine __APPEND2(simulate_check_input,__PREC,__INTSIZE) (tm, s0, &
        tol, prob_low, prob_scale, status)
    !*  SIMULATE_CHECK_INPUT performs input checking for the routines
    !   MARKOV_SIMULATE and MARKOV_SIMULATE_ADVANCED.
    integer, parameter :: INTSIZE = __INTSIZE
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: tm
    integer (INTSIZE), intent(in) :: s0
    real (PREC), intent(in), optional :: tol
        !*  Optional, only used when called from SIMULATE_ADVANCED
    real (PREC), intent(in), optional :: prob_low
        !*  Optional, only used when called from SIMULATE_ADVANCED
    real (PREC), intent(in), optional :: prob_scale
        !*  Optional, only used when called from SIMULATE_ADVANCED
    type (status_t), intent(out) :: status

    integer :: n

    status = NF_STATUS_OK

    n = size(tm,1)

    if (size(tm,1) /= size(tm,2)) goto 100
    if (s0 < 0 .or. s0 > n) goto 100

    call check_nonneg (0.0_PREC, tol, 'tol', status)
    if (status /= NF_STATUS_OK) goto 100

    call check_nonneg (0.0_PREC, prob_low, 'prob_low', status)
    if (status /= NF_STATUS_OK) goto 100

    call check_nonneg (0.0_PREC, prob_scale, 'prob_scale', status)
    if (status /= NF_STATUS_OK) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG
end subroutine



subroutine __APPEND2(simulate,__PREC,__INTSIZE) (tm, s0, seq, status)

    integer, parameter :: INTSIZE = __INTSIZE
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: tm
        !*  Markov process transition matrix
    integer (INTSIZE), intent(in) :: s0
        !*  Initial state for simulation
    integer (INTSIZE), intent(out), dimension(:) :: seq
        !*  Array to store generated sequence of random integers
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:), allocatable :: rnd
    real (PREC), dimension(:,:), allocatable :: cdf
    integer :: i, n, nobs
    integer (INTSIZE) :: ifrom, ito, nstates
    real (PREC) :: x

    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    n = size(tm, 1)
    nobs = size(seq)

    call simulate_check_input (tm, s0, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (nobs == 0) goto 100

    if (n > huge(1_INTSIZE)) then
        ! Number of states implied by trans. matrix size is larger than
        ! what can be represented by integers of kind INTSIZE
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    nstates = int(n, INTSIZE)

    ! draw sequence of random numbers
    allocate (rnd(nobs-1))
    call random_number (rnd)

    ! CDF of z
    allocate (cdf(nstates,nstates))
    cdf(:,1) = tm(:,1)
    do ito = 2, nstates
        cdf(:,ito) = cdf(:,ito-1) + tm(:,ito)
    end do

    seq(1) = s0
    do i = 2, nobs
        ifrom = seq(i-1)
        x = rnd(i-1)

        do ito = 1_INTSIZE, nstates
            if (x <= cdf(ifrom, ito)) exit
        end do

        seq(i) = ito
    end do

100 continue
    if (present(status)) status = lstatus

    if (allocated(rnd)) deallocate(rnd)
    if (allocated(cdf)) deallocate(cdf)

end subroutine



subroutine __APPEND2(simulate_advanced,__PREC,__INTSIZE) (tm, s0, seq, tol, &
            prob_low, prob_scale, tm_sample, status)
    !*  MARKOV_SIMULATE_ADVANCED samples a realization of a Markov process that
    !   satisfies "steady-state" properties even in smaller samples.
    !   Optionally, the routine supports oversampling low-probability
    !   transitions by scaling up the probability of such transitions
    !   (and proportionally scaling down the probabilities of all other
    !   transitions).

    integer, parameter :: INTSIZE = __INTSIZE
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: tm
        !*  Markov process transition matrix
    integer (INTSIZE), intent(in) :: s0
        !*  Initial state for simulation
    integer (INTSIZE), intent(out), dimension(:) :: seq
        !*  Array to store generated sequence of random integers
    real (PREC), intent(in), optional :: tol
        !*  Relative tolerance when comparing sample and population
        !   transition PMFs to determine goodness of fit
        !   (default: 5.0d-2)
    real (PREC), intent(in), optional :: prob_low
        !*  If present, triggers oversampling of those states which have
        !   an unconditional (ergodic) probability below PROB_LOW. Inflows prob.
        !   into such states in given transition matrix are then scaled (up)
        !   by PROB_SCALE, and the transition probabilities to non-oversampled
        !   states are accordingly scaled down.
    real (PREC), intent(in), optional :: prob_scale
        !*  Relative scaling applied to transition probabilities
        !   TM(i,j) <= if I is a state that is to be oversampled and J is
        !   not (default: ignored)
    real (PREC), intent(out), dimension(:,:), optional :: tm_sample
        !*  If present, stores the actual transition matrix used to sample
        !   Markov process (differs from TM only if oversampling is used).
    type (status_t), intent(out), optional :: status

    integer :: nwindow
        ! Window length used to scan for "optimal" subsequences
    integer :: nsim
        ! Total number of time periods to simulate, will be used to pick optimal
        ! subsequences.
    integer :: ifrom, ito
        ! index of first element in current window
    integer :: ifrom_seq
        ! index of current element on seq array
    integer, dimension(:,:), allocatable :: ntrans
        ! transition count in current window
    real (PREC), dimension(:,:), allocatable :: freq_w
        ! Rel. frequency of each transition in sample window
    real (PREC), dimension(:), allocatable :: ergodic_dist
        ! ergodic distribution associated with transition matrix
    real (PREC), dimension(:,:), allocatable :: freq_pop
        ! Rel. frequency of each transition that would be observed in population
    real (PREC), dimension(:,:), allocatable :: tm_s
    real (PREC), dimension(:), allocatable :: pmf
    logical, dimension(:), allocatable :: mask1d
    integer, dimension(:), allocatable :: ios, iall, iother
    integer :: nos, nother
    type (status_t) :: lstatus

    real (PREC) :: eps, lprob_low, sum_low, ltol
    integer (INTSIZE), dimension(:), allocatable :: rint

    integer :: i, n, nobs, k
    logical :: do_oversample

    n = size(tm, 1)
    nobs = size(seq)
    ! No oversampling enforced by default
    eps = 0.0
    lprob_low = 0.0
    do_oversample = .false.
    ltol = 5.0e-2_PREC

    call simulate_check_input (tm, s0, tol, prob_low, prob_scale, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(tol)) ltol = tol

    if (present(prob_low) .and. present(prob_scale)) then
        eps = prob_scale
        lprob_low = prob_low
        do_oversample = .true.
    end if

    ! number of elements in each window (hence nwindow-1 transitions)
    nwindow = nobs - 1
    ! initially simulated sample used to find fitting sub-sequence
    nsim = nwindow * 10

    allocate (ergodic_dist(n))
    allocate (freq_pop(n,n), source=0.0_PREC)

    ! Compute rel. transition frequencies in steady state
    call markov_ergodic_dist (tm, ergodic_dist, inverse=.true.)

    ! compute expected rel. frequency of each transition in population
    forall (i=1:n) freq_pop(i, :) = tm(i, :) * ergodic_dist(i)

    ! Default transition matrix from which to sample
    allocate (tm_s(n,n), source=tm)

    if (do_oversample) then
        ! Adjustment for oversampling: create new "sampling" transition matrix
        ! TM_S which takes into account oversampling of small probabilities.
        allocate (mask1d(n), iall(n))
        call arange (iall)
        mask1d(:) = (ergodic_dist <= lprob_low)

        ! If no element in the ergodic distribution is below the threshold value
        ! don't do anything. Then TM_S and FREQ_POP remain at their initial
        ! values assigned above.
        if (any(mask1d)) then
            allocate (ios(n), iother(n), pmf(n))
            ! Indices of states that need to be oversampled
            nos = count(mask1d)
            ios(:) = 0
            ios(1:nos) = pack(iall, mask1d)
            call setdiff (iall, ios(1:nos), nother, iother)

            do k = 1, nother
                i = iother(k)
                sum_low = sum(tm(i,:), mask=mask1d)
                ! Rescale low-prob. transitions upward, adjust higher-prob.
                ! transitions ownward accordingly.
                where (mask1d)
                    pmf = tm(i,:) * (1.0 + eps)
                else where
                    pmf = tm(i,:) * (1.0 - sum_low * (1.0+eps)) / (1.0 - sum_low)
                end where

                tm_s(i,:) = pmf
            end do

            ! Copy transitions "row" for states with low ergodic prob. as is
            do k = 1, nos
                i = ios(k)
                tm_s(i,:) = tm(i,:)
            end do

            deallocate (ios, iother, pmf)

            ! Recrease ergodic distr. and frequency matrix from updated
            ! transition matrix.
            call markov_ergodic_dist (tm_s, ergodic_dist, inverse=.false.)
            forall (i=1:n) freq_pop(i, :) = tm_s(i, :) * ergodic_dist(i)
        end if

        deallocate (mask1d, iall)
    end if

    ! draw sequence of random numbers
    allocate (rint(nsim))
    call markov_simulate (tm_s, s0, rint, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! ensure that first element is the initial state
    seq(1) = s0
    ! start adding new draws at position 2
    ifrom_seq = 2

    allocate (ntrans(n,n), source=0)
    allocate (freq_w(n,n), source=0.0_PREC)

    ! compute transition count in initial window
10  ifrom = 1
    call count_trans (rint, ifrom, nwindow, ntrans)

    do while (ifrom_seq <= nobs .and. ifrom <= nsim)
        ! compute "goodness of fit" of transition count in window
        freq_w(:,:) = real(ntrans) / (sum(ntrans) + 1.0d0)

        ! also check that the min. prob. transition is at least sampled at
        ! its expected value
        if (accept_sample (freq_pop, freq_w, ltol)) then
            ! add at most k=i-1 elements, ie k-1 transitions, or less if we
            ! cannot store that many
            k = min(nwindow, nobs - ifrom_seq + 1)
            seq(ifrom_seq:ifrom_seq+k-1) = rint(ifrom:ifrom+k-1)
            ! advance initial indices
            ifrom_seq = ifrom_seq + k
            ! allow for last observation to overlap
            ifrom = ifrom + nwindow - 1

            ! not enough elements remaining in simulated pool, exit loop
            if (nsim - ifrom + 1 < nwindow) exit

            ! compute initial transition counts for new non-overlapping window
            call count_trans (rint, ifrom, nwindow, ntrans)
        else
            ! Shift index until we find an initial index that can be concatenated
            ! to whatever is already in seq, ie. has a non-zero trans. prob.
            do while (ifrom < nobs)
                if (tm_s(seq(ifrom_seq-1), rint(ifrom)) > 0) exit
                ifrom = ifrom + 1
            end do
            ! cannot shift window by one element, exit loop
            ito = ifrom + nwindow - 1
            if (ito >= nsim) exit

            ! partially update transition count:
            ! Remove transition at beginning of window
            ntrans(rint(ifrom),rint(ifrom+1)) = ntrans(rint(ifrom),rint(ifrom+1)) - 1
            ! add next transition following at the end of current window
            ! to transition count
            ntrans(rint(ito),rint(ito+1)) = ntrans(rint(ito),rint(ito+1)) + 1
            ! shift window
            ifrom = ifrom + 1
        end if
    end do

    ! check whether we have simulated enough draws, otherwise get a new random
    ! sequence and continue.
    ! Note: this might introduce and infinite loop, add some restriction?
    if (ifrom_seq <= nobs) then
        call markov_simulate (tm_s, s0, rint)
        goto 10
    end if

    if (present(tm_sample)) then
        tm_sample(:,:) = tm_s
    end if

    deallocate (rint)
    deallocate (freq_w, freq_pop, ntrans)
    deallocate (tm_s, ergodic_dist)

100 continue

    if (present(status)) status = lstatus

    contains

    pure subroutine count_trans (x, i0, nmax, num)
        !*  COUNT_TRANS returns the number of observed transitions for each
        !   transition type.
        integer (INTSIZE), intent(in), dimension(:) :: x
        integer, intent(in) :: i0, nmax
        integer, intent(out), dimension(:,:) :: num

        integer :: i

        num = 0
        do i = i0, min(i0 + nmax - 1, size(x) - 1)
            num(x(i), x(i+1)) = num(x(i), x(i+1)) + 1
        end do

    end subroutine

    pure function accept_sample (x, xhat, tol) result(res)
        !*  ACCEPT_SAMPLE determines whether a simulated sample satisfies
        !   the acceptance criteria.
        real (PREC), intent(in), dimension(:,:) :: x
            !*  Matrix containing the steady-state transition PMF
        real (PREC), intent(in), dimension(:,:) :: xhat
            !*  Matrix containing sample transition PMF
        real (PREC), intent(in) :: tol
            !*  Acceptable error, relative to steady-state PMF
        logical :: res

        real (PREC) :: maxdiff
        real (PREC), dimension(:,:), allocatable :: rdiff

        allocate (rdiff(size(x,1),size(x,2)))

        where (x > 0)
            rdiff = abs(xhat - x) / x
        else where
            rdiff = 0.0d0
        end where
        maxdiff = maxval(rdiff)

        deallocate (rdiff)

        res = (maxdiff <= tol)
    end function

end subroutine

