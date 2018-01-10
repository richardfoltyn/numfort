
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

    if (nobs == 0) goto 100
    if (size(tm,2) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

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
