real (PREC), intent(out), dimension(:) :: x
real (PREC), intent(in) :: xfrom, xto

real (PREC) :: step
integer :: i, n

n = size(x)
step = (xto - xfrom) / (n - 1.0_PREC)

x(1) = xfrom
do i = 2, n-1
    x(i) =  xfrom + (i-1) * step
end do

! avoid rounding errors in end point
x(n) = xto
