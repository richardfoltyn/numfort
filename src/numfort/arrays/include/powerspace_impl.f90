real (PREC), intent(out), dimension(:) :: x
    !*  Array to store "power-spaced" sequence of points
real (PREC), intent(in) :: xmin
    !*  Starting value
real (PREC), intent(in) :: xmax
    !*  Endpoint value
real (PREC), intent(in) :: pow
    !*  Exponent used to create power-spaced sequence

integer :: n, i
real (PREC) :: slope

n = size(x)
call linspace (x, 0.0_PREC, 1.0_PREC)

slope = xmax - xmin
do i = 1, n
    x(i) = xmin + slope * x(i) ** pow
end do

! Explicitly set boundary values to eliminate any rounding issues
x(1) = xmin
x(n) = xmax
