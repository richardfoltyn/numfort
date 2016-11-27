! Implementation for 1D linear interpolation

real (PREC), intent(in), dimension(:) :: x
    !!  The x-coordinates of the interpolated values.
real (PREC), intent(in), dimension(:) :: xp
    !!  The x-coordinates of data points in increasing order.
real (PREC), intent(in), dimension(:) :: fp
    !!  The y-coordinates of data points, same length as xp.
logical, intent(in), optional :: ext
    !!  If present and true, function values for x-coordinates outside of
    !!  boundary values of xp are linearly extrapolated. Default is `.false.`.
real (PREC), intent(in), optional :: left
    !!  Value to return if x < xp(1) and extrapolation is disabled (`ext=.false.`).
    !!  Default is fp(1)
real (PREC), intent(in), optional :: right
    !!  Value to return if x > xp(-1) and extrapolation disabled (`ext=.false.`).
    !!  Default is fp(-1)
real (PREC), intent(out), dimension(:) :: fx
    !!  Array of interpolated values, same shape as x.

integer (INTSIZE) :: n, i
logical :: lext
real (PREC) :: lright, lleft

n = size (x)
lext = .false.

if (present(ext)) lext = ext
if (present(left)) then
    lleft = left
else
    lleft = fp(1)
end if

if (present(right)) then
    lright = right
else
    lright = fp(size(fp))
end if

do i = 1, n
    ! call scalar implementation
    call interp_linear_impl (x(i), xp, fp, fx(i), lext, lleft, lright)
end do
