!*  Implementation code of DELETE routine for 2d-array that is indepenent
!   of argument type/kind.

    integer, intent(in), dimension(:), target :: idx
        !*  Index of elements which should be deleted from ARR.
    integer, intent(in) :: dim
        !*  The dimension along with elements should be deleted.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted by
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT, dimension DIM,
        !   that holds valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer :: narr, nout, nidx, i, j, ioffset, k, m, ioffset_a, nkeep
    integer, dimension(:), pointer :: uidx
    integer, dimension(:), allocatable :: ikeep, iall
    logical :: lsorted
    type (status_t) :: lstatus

    lstatus = NF_STATUS_INVALID_ARG
    ! Offset on OUT array
    ioffset = 0

    ! Ignore the DIM argument in the 1d-implementation, but exit if it has
    ! an invalid value.
    if (dim < 1 .or. dim > 2) goto 100

    lsorted = .false.
    if (present(sorted)) lsorted = sorted

    nidx = size(idx)
    ! If required, sort index array and retain only unique elements
    if (.not. lsorted) then
        allocate (uidx(nidx))
        call unique (idx, nidx, uidx)
    else
        uidx => idx(1:nidx)
    end if

    narr = size(arr, dim)
    nout = size(out, dim)
    ! Size of "other" dimension on input array
    m = size(arr, 3-dim)

    ! Check whether there is enough space along the dimension where deletion occurs
    if (nout < (narr - nidx)) goto 100
    ! Check space in "other" dimension
    if (size(out,3-dim) < m) goto 100

    ! Handle empty index array separately
    if (nidx == 0) then
        m = size(arr, 1)
        ! DIM is irrelevant, just copy over data column-wise
        forall (i=1:size(arr,2)) out(1:m,i) = arr(1:m,i)
        lstatus = NF_STATUS_OK
        ioffset = narr
        goto 100
    end if

    ! Check that indices are valid
    if (uidx(1) < 1) goto 100
    ! Check whether any index exceeds upper bound before beginning to copy values.
    ! Sufficient to check last index.
    if (uidx(nidx) > narr) goto 100

    if (dim == 1) then
        ! DIM = 1, so selected rows should be deleted.
        nkeep = narr - nidx
        ! Find indices or all rows that should be kept so we don't have to
        ! repeat this for every column.
        allocate (ikeep(nkeep), iall(narr))
        call arange (iall)
        call setdiff (iall, uidx(:nidx), nkeep, ikeep, assume_unique=.true.)

        ! Copy by column; within each column, copy over elements that should
        ! not be deleted.
        do j = 1, size(arr, 2)
            do i = 1, nkeep
                k = ikeep(i)
                out(i,j) = arr(k,j)
            end do
        end do

        ioffset = nkeep
    else
        ! DIM = 2, so colums should be deleted. In this case we just copy
        ! over column by column, omitting those corresponding to values in
        ! UIDX.

        ! Copy first segment in ARR that is located before the first index
        ! to be dropped.
        ioffset_a = 0
        k = uidx(1)-1
        forall (i=1:k) out(1:m,ioffset+i) = arr(1:m,ioffset_a+i)
        ioffset = k

        do i = 1, nidx-1

            ! Not enough space in OUT to store any additinal data
            if (ioffset >= nout) goto 10

            ioffset_a = uidx(i)
            k = uidx(i+1) - ioffset_a - 1
            forall (j=1:k) out(1:m,ioffset+j) = arr(1:m,ioffset_a+j)
            ioffset = ioffset + k
        end do

        ! Copy last segment that with element indices uidx(nidx),...,narr
        ioffset_a = uidx(nidx)
        k = narr - ioffset_a
        forall (i=1:k) out(1:m,ioffset+i) = arr(1:m,ioffset_a+i)
        ioffset = ioffset + k
    end if

10  continue

    ! At this point all elements from ARR that should be retained were copied
    ! over to OUT to the extend possible.
    lstatus = NF_STATUS_OK

100 continue
    if (present(status)) status = lstatus
    if (present(n)) n = ioffset

    if (.not. lsorted .and. associated(uidx)) then
        deallocate (uidx)
    end if

    if (allocated(ikeep)) deallocate (ikeep)
    if (allocated(iall)) deallocate (iall)
