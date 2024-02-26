

subroutine simulate_check_input (tm, s0, rtol, atol, prob_low, prob_scale, transposed, &
        status)
    !*  SIMULATE_CHECK_INPUT performs input checking for the routines
    !   MARKOV_SIMULATE and MARKOV_SIMULATE_ADVANCED.
    real (PREC), intent(in), dimension(:,:) :: tm
    integer (INTSIZE), intent(in) :: s0
    real (PREC), intent(in), optional :: rtol
        !*  Optional, only used when called from SIMULATE_ADVANCED
    real (PREC), intent(in), optional :: atol
        !*  Optional, only used when called from SIMULATE_ADVANCED
    real (PREC), intent(in), optional :: prob_low
        !*  Optional, only used when called from SIMULATE_ADVANCED
    real (PREC), intent(in), optional :: prob_scale
        !*  Optional, only used when called from SIMULATE_ADVANCED
    logical, intent(in), optional :: transposed
    type (status_t), intent(out) :: status

    integer :: n

    status = NF_STATUS_OK

    n = size(tm,1)

    if (size(tm,1) /= size(tm,2)) goto 100
    if (s0 < 1 .or. s0 > n) goto 100

    ! Check for valid transition matrix
    if (.not. is_trans_matrix (tm, transposed=transposed)) goto 100

    call check_nonneg (0.0_PREC, rtol, 'rtol', status)
    if (status /= NF_STATUS_OK) goto 100

    call check_nonneg (0.0_PREC, atol, 'atol', status)
    if (status /= NF_STATUS_OK) goto 100

    call check_range (0.0_PREC, prob_low, 0.0_PREC, 1.0_PREC, 'prob_low', status)
    if (status /= NF_STATUS_OK) goto 100

    call check_nonneg (0.0_PREC, prob_scale, 'prob_scale', status)
    if (status /= NF_STATUS_OK) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG
end subroutine



subroutine simulate (transm, s0, seq, transposed, status)
    real (PREC), intent(in), dimension(:,:), contiguous, target :: transm
        !*  Markov process transition matrix with each element (i,j)
        !   representing the transition probability Prob[x'=j | x=i].
    integer (INTSIZE), intent(in) :: s0
        !*  Initial state for simulation
    integer (INTSIZE), intent(out), dimension(:), contiguous :: seq
        !*  Array to store generated sequence of random integers
    logical, intent(in), optional :: transposed
        !*  If present and true, TM is assumed to be in transposed format,
        !   ie. each element (i,j) represents the probability Prob[x'=i | x = j]
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:), allocatable :: rnd
    real (PREC), dimension(:,:), allocatable :: cdf
    integer :: n, nobs, iobs
    integer (INTSIZE) :: ifrom, ito, nstates, i
    real (PREC) :: x
    logical :: ltransposed

    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    n = size(transm, 1)
    nobs = size(seq)
    ltransposed = .false.
    if (present(transposed)) ltransposed = transposed

    call simulate_check_input (transm, s0, transposed=transposed, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (nobs == 0) goto 100

    if (n > huge(1_INTSIZE)) then
        ! Number of states implied by trans. matrix size is larger than
        ! what can be represented by integers of kind INTSIZE
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    nstates = int(n, INTSIZE)

    ! Create transposed matrix if applicable as that one can be traversed
    ! faster in the loop below
    if (ltransposed) then
        allocate (cdf(nstates,nstates), source=transm)
    else
        allocate (cdf(nstates,nstates), source=transpose(transm))
    end if

    ! draw sequence of random numbers
    allocate (rnd(nobs-1))
    call random_number (rnd)

    ! CDF of z
    ! At this points CDF contains the (transposed) transition matrix
    do i = 1, nstates
        do ito = 2, nstates
           cdf(ito,i) = cdf(ito-1,i) + cdf(ito,i)
        end do
    end do

    seq(1) = s0
    do iobs = 2, nobs
        ifrom = seq(iobs-1)
        x = rnd(iobs-1)

        do ito = 1_INTSIZE, nstates
            if (x <= cdf(ito,ifrom)) exit
        end do

        seq(iobs) = ito
    end do

100 continue

    if (present(status)) status = lstatus

    if (allocated(rnd)) deallocate(rnd)
    if (allocated(cdf)) deallocate(cdf)

end subroutine



subroutine simulate_advanced (tm, s0, seq, atol, rtol, prob_low, prob_scale, &
        tm_sample, transposed, status)
    !*  MARKOV_SIMULATE_ADVANCED samples a realization of a Markov process that
    !   satisfies "steady-state" properties even in smaller samples.
    !   Optionally, the routine supports oversampling low-probability
    !   transitions by scaling up the probability of such transitions
    !   (and proportionally scaling down the probabilities of all other
    !   transitions).
    real (PREC), intent(in), dimension(:,:), contiguous, target :: tm
        !*  Markov process transition matrix
    integer (INTSIZE), intent(in) :: s0
        !*  Initial state for simulation
    integer (INTSIZE), intent(out), dimension(:), contiguous :: seq
        !*  Array to store generated sequence of random integers
    real (PREC), intent(in), optional :: rtol
        !*  Relative tolerance when comparing sample and population
        !   transition PMFs to determine goodness of fit
        !   (default: 5.0e-2)
    real (PREC), intent(in), optional :: atol
        !*  Absolute tolerance when comparing sample and population transition
        !   PMFs to determine goodness of fit
        !   (default: 1.0e-3)
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
    logical, intent(in), optional :: transposed
        !*  If present and true, TM is assumed to be in transposed format,
        !   ie. each element (i,j) represents the probability Prob[x'=i | x = j]
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
    real (PREC), dimension(:,:), allocatable :: ltm
    type (status_t) :: lstatus

    real (PREC) :: eps, lprob_low, lrtol, latol
    integer (INTSIZE), dimension(:), allocatable :: rint

    integer :: i, n, nobs, k
    logical :: do_oversample, ltransposed

    n = size(tm, 1)
    nobs = size(seq)
    ! No oversampling enforced by default
    eps = 0.0
    lprob_low = 0.0
    do_oversample = .false.
    lrtol = 5.0e-2_PREC
    latol = 1.0e-3_PREC

    call simulate_check_input (tm, s0, rtol, atol, prob_low, prob_scale, &
        transposed, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ltransposed = .false.
    if (present(transposed)) ltransposed = transposed
    if (present(rtol)) lrtol = rtol
    if (present(atol)) latol = atol

    if (present(prob_low) .and. present(prob_scale)) then
        eps = prob_scale
        lprob_low = prob_low
        do_oversample = .true.
    end if

    ! Create tranposed transition matrix
    if (ltransposed) then
        allocate (ltm(n,n), source=tm)
    else
        allocate (ltm(n,n), source=transpose(tm))
    end if

    ! number of elements in each window (hence nwindow-1 transitions)
    nwindow = nobs - 1
    ! initially simulated sample used to find fitting sub-sequence
    nsim = nwindow * 10

    allocate (ergodic_dist(n))
    allocate (freq_pop(n,n), source=0.0_PREC)

    if (do_oversample) then
        ! Adjust transition matrix in place for oversampling
        call oversample_tm (ltm, lprob_low, prob_scale)
    end if

    ! ergodic distr. and transition frequency matrix
    call markov_ergodic_dist (ltm, ergodic_dist, transposed=.true., inverse=.false.)
    do i = 1, n
        freq_pop(:,i) = ltm(:,i) * ergodic_dist(i)
    end do

    ! draw sequence of random numbers
    allocate (rint(nsim))
    call markov_simulate (ltm, s0, rint, transposed=.true., status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! ensure that first element is the initial state
    seq(1) = s0
    ! start adding new draws at position 2
    ifrom_seq = 2

    allocate (ntrans(n,n), source=0)
    allocate (freq_w(n,n), source=0.0_PREC)

    ! compute transition count in initial window
10  ifrom = 1
    call count_trans (rint(ifrom:ifrom+nwindow), ntrans)

    do while (ifrom_seq <= nobs .and. ifrom <= nsim)
        ! compute "goodness of fit" of transition count in window
        freq_w(:,:) = real(ntrans, PREC) / (sum(ntrans) + 1.0_PREC)

        ! also check that the min. prob. transition is at least sampled at
        ! its expected value
        if (all_close_fast (freq_w, freq_pop, rtol=lrtol, atol=latol)) then
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
            call count_trans (rint(ifrom:nwindow), ntrans)
        else
            ! Shift index until we find an initial index that can be concatenated
            ! to whatever is already in seq, ie. has a non-zero trans. prob.
            do while (ifrom < nobs)
                if (ltm(rint(ifrom), seq(ifrom_seq-1)) > 0) exit
                ifrom = ifrom + 1
            end do
            ! cannot shift window by one element, exit loop
            ito = ifrom + nwindow - 1
            if (ito >= nsim) exit

            ! partially update transition count:
            ! Remove transition at beginning of window
            ntrans(rint(ifrom+1),rint(ifrom)) = ntrans(rint(ifrom+1),rint(ifrom)) - 1
            ! add next transition following at the end of current window
            ! to transition count
            ntrans(rint(ito+1),rint(ito)) = ntrans(rint(ito+1),rint(ito)) + 1
            ! shift window
            ifrom = ifrom + 1
        end if
    end do

    ! check whether we have simulated enough draws, otherwise get a new random
    ! sequence and continue.
    ! Note: this might introduce and infinite loop, add some restriction?
    if (ifrom_seq <= nobs) then
        call markov_simulate (ltm, s0, rint, transposed=.true.)
        goto 10
    end if

    ! Write actual sampling transition matrix to user-provided output array
    if (present(tm_sample)) then
        if (ltransposed) then
            tm_sample(:,:) = ltm
        else
            tm_sample(:,:) = transpose(ltm)
        end if
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine oversample_tm (transm, prob_low, prob_scale)
    real (PREC), intent(inout), dimension(:,:), contiguous :: transm
        !*  Transposed input transition matrix. Contains oversampled
        !   transposed transition matrix on exit.
    real (PREC), intent(in) :: prob_low
    real (PREC), intent(in) :: prob_scale

    logical, dimension(:), allocatable :: mask1d
    integer, dimension(:), allocatable :: ios, iall, iother
    real (PREC), dimension(:), allocatable :: ergodic_dist, pmf
    real (PREC), dimension(:,:), allocatable :: transm_s
    integer :: n, nos, nother
    integer :: i, k
    real (PREC) :: sum_low

    n = size(transm, 1)

    allocate (ergodic_dist(n))

    ! Compute rel. transition frequencies in steady state
    call markov_ergodic_dist (transm, ergodic_dist, transposed=.true., inverse=.true.)

    ! Adjustment for oversampling: create new "sampling" transition matrix
    ! TM_S which takes into account oversampling of small probabilities.
    allocate (mask1d(n), iall(n))
    call arange (iall)
    mask1d(:) = (ergodic_dist <= prob_low)

    ! If no element in the ergodic distribution is below the threshold value
    ! don't do anything.
    if (any(mask1d)) then
        allocate (transm_s(n,n), source=0.0_PREC)
        allocate (ios(n), iother(n), pmf(n))
        ! Indices of states that need to be oversampled
        nos = count(mask1d)
        ios(:) = 0
        ios(1:nos) = pack(iall, mask1d)
        call setdiff (iall, ios(1:nos), nother, iother)

        do k = 1, nother
            i = iother(k)
            sum_low = sum(transm(:,i), mask=mask1d)
            ! Rescale low-prob. transitions upward, adjust higher-prob.
            ! transitions ownward accordingly.
            where (mask1d)
                pmf = transm(:,i) * (1.0 + prob_scale)
            else where
                pmf = transm(:,i) * (1.0 - sum_low * (1.0_PREC+prob_scale)) / (1.0 - sum_low)
            end where

            transm_s(:,i) = pmf
        end do

        ! Copy transitions "row" for states with low ergodic prob. as is
        do k = 1, nos
            i = ios(k)
            transm_s(:,i) = transm(:,i)
        end do

        ! Write back oversampled transition matrix
        transm(:,:) = transm_s(:,:)
    end if

end subroutine



pure subroutine count_trans (x, num)
    !*  COUNT_TRANS returns the number of observed transitions for each
    !   transition type.
    integer (INTSIZE), intent(in), dimension(:), contiguous :: x
    integer, intent(out), dimension(:,:), contiguous :: num
        !*  Matrix containing the number of transitions for each possible
        !   combination of origin and destrination values, in transposed format.
        !   Each element (i,j) thus contains the number of transitions
        !   for j to i.

    integer :: i, n
    integer (INTSIZE) :: ifrom,  ito

    n = size(x)
    num = 0

    do i = 1, n - 1
        ifrom = x(i)
        ito = x(i+1)
        num(ito, ifrom) = num(ito, ifrom) + 1
    end do

end subroutine

