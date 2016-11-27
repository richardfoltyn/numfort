! Wrapper around interp_linear for scalar x, fx

real (PREC), intent(in) :: x
    !!  The x-coordinate of the interpolated value.
real (PREC), intent(in), dimension(:) :: xp
    !!  The x-coordinates of data points in increasing order.
real (PREC), intent(in), dimension(:) :: fp
    !!  The y-coordinates of data points, same length as xp.
logical, intent(in), optional :: ext
    !!  If present and true, function values for x-coordinates outside of
    !!  boundary values of xp are linearly extrapolated. Default is `.false.`.
real (PREC), intent(in), optional :: left
    !!  Value to return if x < xp(1) and extrapolation is disabled (`ext=.false.`).
    !!  Default is fp(1).
real (PREC), intent(in), optional :: right
    !!  Value to return if x > xp(-1) and extrapolation disabled (`ext=.false.`).
    !!  Default is fp(-1).
real (PREC), intent(out) :: fx
    !!  Interpolated value

logical :: lext
real (PREC) :: lright, lleft

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

call interp_linear_impl (x, xp, fp, fx, lext, lleft, lright)
