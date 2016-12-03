! Adapt source code from RANDOM library for Fortran 2008

integer (INTSIZE), intent(out), dimension(:) :: x
    !!  Array where permutation should be stored.
integer (INTSIZE), intent(in), optional :: low
    !!  If present, create permutation of set {low, ..., N} where
    !!  N = low + size(x) - 1 (default: 1)

integer (INTSIZE) :: i, j, k, n, offset
real (real64) :: wk

if (present(low)) then
    offset = low - 1
else
    offset = 0
end if

n = size(x)
forall (i=1:n) x(i) = i + offset

! Starting at the end, swap the current last indicator with one
! randomly chosen from those preceeding it.
do i = n, 2, -1
    call random_number (wk)
    j = 1 + i * wk
    if (j < i) then
        k = x(i)
        x(i) = x(j)
        x(j) = k
    end if
end do
