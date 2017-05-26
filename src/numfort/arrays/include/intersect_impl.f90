integer (INTSIZE), intent(out) :: n
    !*  On exit, contains the number of elements in intersection
logical, intent(in), optional :: assume_unique
    !*  If present and true, arrays arr1 and arr2 are assumed to contain no
    !   duplicates. Elements need not be sorted.

logical :: lassume_unique

integer :: m1, m2, nidx1, nidx2, nwork, i
integer, dimension(:), allocatable :: idx1, idx2, iwork


lassume_unique = .false.
if (present(assume_unique)) lassume_unique = assume_unique

m1 = size(arr1)
m2 = size(arr2)


if (lassume_unique) then
    nwork = m1 + m2
    allocate (work(nwork), iwork(nwork))
    forall (i=1:m1) work(i) = arr1(i)
    forall (i=1:m2) work(m1+i) = arr2(i)
else
    ! Need to determine unique elements first
    allocate (idx1(m1), idx2(m2))
    call unique (arr1, nidx1, idx=idx1)
    call unique (arr2, nidx2, idx=idx2)

    nwork = nidx1 + nidx2

    ! Pool unique elements in work array
    allocate (work(nwork), iwork(nwork))
    forall (i=1:nidx1) work(i) = arr1(idx1(i))
    forall (i=1:nidx2) work(nidx1+i) = arr2(idx2(i))

    deallocate (idx1, idx2)
end if

! Sort pooled arrays of unique elements such that elements that are present
! in both arrays will end right next to each other.
call argsort (work, iwork)

n = 0
i = 1
do while (i < nwork)
    if (work(iwork(i)) == work(iwork(i+1))) then
        n = n + 1
        if (n <= size(res)) res(n) = work(iwork(i))
        ! If element of intersection found, we can skip forward one additional
        ! element, as we are guaranteed that elements at i+1 and i+2 will
        ! be different.
        i = i + 1
    end if

    i = i + 1
end do

deallocate (work, iwork)
