real (PREC), intent(out), dimension(:) :: x
real (PREC), intent(in) :: xfrom, xto

real (PREC) :: step
integer :: i, n

n = size(x)
step = (xto - xfrom) / (n - 1)

x(1) = xfrom
do i = 1, n-2
    x(i) =  xfrom + i * step
end do

! avoid rounding errors in end point
x(n) = xto
