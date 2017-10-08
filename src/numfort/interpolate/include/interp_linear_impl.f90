!*  Implements linear interpolation for scalar input/return value
!   and no optional arguments.
!   Should be called from wrapper routines doing the optional argument handling,
!   etc.

real (PREC), intent(in) :: x
    !*  x-coordinate of the interpolated value
real (PREC), intent(in), dimension(:) :: xp
    !*  x-coordinates of data points in increasing order
real (PREC), intent(in), dimension(:) :: fp
    !*  y-coordinates of the data points, same length as xp.
integer (NF_ENUM_KIND), intent(in) :: ext
    !*  Defines behavious in case function should be evaluated at point x
    !   outside of the range specified by [xp(1), xp(size(xp))].
    !   Adminissible constants are 
    !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
    !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
    !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
    !           are used, respectively.
    !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and 
    !           fp(size(xp)) for all x > xp(size(xp)).
    !       - NF_INTERP_EVAL_ZERO: Returns zero for all x outside of
    !           admissible range.
    !       - For all other values, the function value at x is extrapolated.
real (PREC), intent(in) :: left
    !*  Value to return if x < xp(1)
real (PREC), intent(in) :: right
    !*  Value to return if x > xp(size(xp))
real (PREC), intent(out) :: fx
    !*  Interpolated value

integer (INTSIZE) :: ilb, iub
real (PREC) :: slope
integer :: np

np = size(xp)

if (x < xp(1)) then
    select case (ext)
    case (NF_INTERP_EVAL_CONST)
        fx = left
        return
    case (NF_INTERP_EVAL_BOUNDARY)
        fx = fp(1)
        return
    case (NF_INTERP_EVAL_ZERO)
        fx = 0.0_PREC
        return
    end select
else if (x > xp(np)) then
    select case (ext)
    case (NF_INTERP_EVAL_CONST)
        fx = right
        return
    case (NF_INTERP_EVAL_BOUNDARY)
        fx = fp(np)
        return
    case (NF_INTERP_EVAL_ZERO)
        fx = 0.0_PREC
        return
    end select
end if

! default: At this point there is either a bracketing interval, or we extrapolate
! values outside of domain
ilb = interp_find (x, xp)
iub = ilb + 1
slope = (x - xp(ilb)) / (xp(iub) - xp(ilb))
fx = (1.0_PREC-slope) * fp(ilb) + slope * fp(iub)
