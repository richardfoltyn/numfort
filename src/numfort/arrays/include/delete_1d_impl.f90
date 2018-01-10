!*  Implementation code of DELETE routine for 1d-array that is indepenent
!   of argument type/kind.

    integer, intent(in), dimension(:), target :: idx
        !*  Index of elements which should be deleted from ARR.
    integer, intent(in), optional :: dim
        !*  The dimension along with elements should be deleted. Ignored for
        !   1d-arrays. An error is raised if this argument is present and DIM /= 1.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted by
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT that holds
        !   valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer :: narr, nout, nidx, i, ifrom, ito, ioffset, k
    integer, dimension(:), pointer :: uidx
    logical :: lsorted
    type (status_t) :: lstatus

    nullify (uidx)

    lsorted = .false.
    if (present(sorted)) lsorted = sorted

    lstatus = NF_STATUS_INVALID_ARG
    ! Offset on OUT array
    ioffset = 0

    ! Ignore the DIM argument in the 1d-implementation, but exit if it has
    ! an invalid value.
    if (present(dim)) then
        if (dim /= 1) goto 100
    end if

    narr = size(arr)
    nout = size(out)
    nidx = size(idx)

    ! Handle empty index array separately
    if (nidx == 0) then
        ifrom = 1
        ito = min(narr, nout)
        out(ifrom:ito) = arr(ifrom:ito)
        lstatus = NF_STATUS_OK
        ioffset = ito
        goto 100
    end if

    ! If required, sort index array and retain only unique elements
    if (.not. lsorted) then
        allocate (uidx(nidx))
        call unique (idx, nidx, uidx)
    else
        uidx => idx(1:nidx)
    end if

    ! Check that indices are valid
    if (uidx(1) < 1) goto 100
    ! Check whether any index exceeds upper bound before beginning to copy values.
    ! Sufficient to check last index.
    if (uidx(nidx) > narr) goto 100

    ! Copy first segment in ARR that is located before the first index
    ! to be dropped.
    ifrom = 1
    k = min(uidx(1)-1, nout)
    ito = ifrom + k - 1
    out(ioffset+1:ioffset+k) = arr(ifrom:ito)
    ioffset = k

    do i = 1, nidx-1

        ! Not enough space in OUT to store any additinal data
        if (ioffset >= nout) goto 10

        ifrom = uidx(i) + 1
        k = min(uidx(i+1)-ifrom, nout-ioffset)
        ito = ifrom + k - 1

        out(ioffset+1:ioffset+k) = arr(ifrom:ito)
        ioffset = ioffset + k
    end do

    ! Copy last segment that with element indices uidx(nidx)+1,...,narr
    ifrom = uidx(nidx) + 1
    k = min(narr-ifrom+1, nout-ioffset)
    ito = ifrom + k - 1
    out(ioffset+1:ioffset+k) = arr(ifrom:ito)
    ioffset = ioffset + k

10  continue

    ! At this point all elements from ARR that should be retained were copied
    ! over to OUT
    lstatus = NF_STATUS_OK

100 continue
    if (present(status)) status = lstatus
    if (present(n)) n = ioffset

    if (.not. lsorted .and. associated(uidx)) then
        deallocate (uidx)
    end if
