integer (INTSIZE), intent(out) :: n
    !*  On exit, contains the number of elements in the set C = A\B
integer (INTSIZE), intent(out), dimension(:), optional :: idx
    !*  On exit, contains those indices from arr1 which are in the
    !   set difference C = A\B. If idx is not large enough not hold all
    !   elements of C, only the first size(idx) elements are stored and
    !   n > size(idx).
logical, intent(in), optional :: assume_unique
    !*  If present and true, arrays arr1 and arr2 are assumed to contain
    !   no duplicates. Elements need not be sorted.

logical :: lassume_unique
integer :: i, j, k, offset1, offset2, ito
integer :: m1, m2, nidx1, nidx2, ndiff_max, nidx_max
integer, dimension(:), allocatable :: idx1, idx2

lassume_unique = .false.
if (present(assume_unique)) lassume_unique = assume_unique

m1 = size(arr1)
m2 = size(arr2)

allocate (idx1(m1), idx2(m2))

if (lassume_unique) then
    call argsort (arr1, idx1)
    call argsort (arr2, idx2)
    nidx1 = m1
    nidx2 = m2
else
    call unique (arr1, nidx1, idx=idx1)
    call unique (arr2, nidx2, idx=idx2)
end if

ndiff_max = 0
nidx_max = 0
if (present(diff)) ndiff_max = size(diff)
if (present(idx)) nidx_max = size(idx)

if (nidx1 > 0 .and. nidx2 > 0) then
    ! Number of elements in A\B found so far
    n = 0
    offset1 = 1

    ! Discard all elements in B that are smaller than the smallest
    ! value in A.
    do offset2 = 1, nidx2
        if (arr2(idx2(offset2)) >= arr1(idx1(1))) exit
    end do

    if (offset2 <= nidx2) then
        do i = offset2, nidx2
            val = arr2(idx2(i))
            do j = offset1, nidx1
                k = idx1(j)
                if (arr1(k) == val) then
                    offset1 = j + 1
                    exit
                end if

                n = n + 1
                if (n <= ndiff_max) diff(n) = arr1(k)
                if (n <= nidx_max) idx(n) = k
            end do
        end do
    end if

    ! Append all remaining elements in arr1 that
    ! offset1 now points to element right before the first
    ! of the remaining elements that need to be inserted into result.
    offset1 = offset1 - 1
    ! Need to insert nidx1-offset1 elements, but limit by capacity
    ! constraint.
    ito = min(ndiff_max - n, nidx1 - offset1)
    ! n contains number of elements inserted so far.
    forall (i=1:ito) diff(n+i) = arr1(idx1(offset1+i))

    ito = min(nidx_max - n, nidx1 - offset1)
    forall (i=1:ito) idx(n+i) = idx1(offset1+i)

    n = n + (nidx1 - offset1)

else
    ! Set B is empty, just copy over sorted, unique set A to the extent possible
    n = nidx1
    forall (i=1:min(n,ndiff_max)) diff(i) = arr1(idx1(i))
    forall (i=1:min(n,nidx_max)) idx(i) = idx1(i)
end if

deallocate (idx1, idx2)
