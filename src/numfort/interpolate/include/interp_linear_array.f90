! Implementation for 1D linear interpolation

real (PREC), intent(in), dimension(:) :: x
    !*  The x-coordinates of the interpolated values.
real (PREC), intent(in), dimension(:) :: xp
    !*  The x-coordinates of data points in increasing order.
real (PREC), intent(in), dimension(:) :: fp
    !*  The y-coordinates of data points, same length as xp.
integer (NF_ENUM_KIND), intent(in), optional :: ext
    !*  Defines behavious in case function should be evaluated at point x
    !   outside of the range specified by [xp(1), xp(size(xp))].
    !   Adminissible constants are 
    !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
    !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
    !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and 
    !           fp(size(xp)) for all x > xp(size(xp)).
    !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
    !           are used, respectively.
    !       - NF_INTERP_EVAL_ZERO: Returns zero for all x outside of
    !           admissible range.
    !       - For all other values, the function value at x is extrapolated.
real (PREC), intent(in), optional :: left
    !*  Value to return if x < xp(1) and extrapolation is disabled (`ext=.false.`).
    !   Default is fp(1)
real (PREC), intent(in), optional :: right
    !*  Value to return if x > xp(-1) and extrapolation disabled (`ext=.false.`).
    !   Default is fp(-1)
real (PREC), intent(out), dimension(:) :: fx
    !*  Array of interpolated values, same shape as x.

integer (INTSIZE) :: n, i
integer (NF_ENUM_KIND) :: lext
real (PREC) :: lright, lleft

n = size (x)
lext = NF_INTERP_EVAL_EXTRAPOLATE

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
