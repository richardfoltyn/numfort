! Adapt source code from RANDOM library for Fortran 2008

integer (int64) :: i, j, n
real (real64) :: wk

n = size (x)

! Starting at the end, swap the current last indicator with one
! randomly chosen from those preceeding it.
do i = n, 2, -1
    call random_number (wk)
    j = 1 + int(i * wk, int64)
    if (j < i) then
        xi = x(i)
        x(i) = x(j)
        x(j) = xi
    end if
end do
