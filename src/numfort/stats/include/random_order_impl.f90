! Adapt source code from RANDOM library for Fortran 2008

integer (INTSIZE), intent(out), dimension(:) :: x

integer (INTSIZE) :: i, j, k, n
real (real64) :: wk

n = size(x)
forall (i=1:n) x(i) = i

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
