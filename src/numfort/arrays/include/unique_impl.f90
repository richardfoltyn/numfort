

integer (INTSIZE), intent(out) :: n
integer (INTSIZE), intent(in out), dimension(:), optional :: idx

integer (INTSIZE), dimension(:), allocatable :: lidx
integer :: m, k, i

m = size(arr)
allocate (lidx(m))

call unirnk (arr, lidx, n)

if (present(uarr)) then
    k = min(n, size(uarr))
    do i = 1, k
        uarr(i) = arr(lidx(i))
    end do
end if

if (present(idx)) then
    k = min(n, size(idx))
    idx(1:k) = lidx(1:k)
    idx(k+1:) = 0
end if

deallocate (lidx)
