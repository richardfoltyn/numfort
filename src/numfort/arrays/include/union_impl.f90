
integer (INTSIZE), intent(out) :: n
    !*  On exit, contains the number of elements in intersection
logical, intent(in), optional :: assume_unique
    !*  If present and true, arrays arr1 and arr2 are assumed to contain no
    !   duplicates. Elements need not be sorted.

logical :: lassume_unique

integer :: m1, m2, nidx1, nidx2, nwork, i
integer, dimension(:), allocatable :: idx1, idx2


lassume_unique = .false.
if (present(assume_unique)) lassume_unique = assume_unique

m1 = size(arr1)
m2 = size(arr2)


if (lassume_unique) then
    nwork = m1 + m2
    allocate (work(nwork))
    do i = 1, m1
        work(i) = arr1(i)
    end do
    do i = 1, m2
        work(m1+i) = arr2(i)
    end do
else
    ! Need to determine unique elements first
    allocate (idx1(m1), idx2(m2))
    call unique (arr1, nidx1, idx=idx1)
    call unique (arr2, nidx2, idx=idx2)

    nwork = nidx1 + nidx2

    ! Pool unique elements in work array
    allocate (work(nwork))
    do i = 1, nidx1
        work(i) = arr1(idx1(i))
    end do
    do i = 1, nidx2
        work(nidx1+i) = arr2(idx2(i))
    end do

    deallocate (idx1, idx2)
end if

! Sort pooled arrays of unique elements such that elements that are present
! in both arrays will end right next to each other.
call unique (work, n, uarr=res)


deallocate (work)
