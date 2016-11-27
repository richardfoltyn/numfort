! Implementation of linear interpolation for scalar input/return value
! and no optional arguments.
! Should be called from wrapper routines doing the optional argument handling,
! etc.

!>  x-coordinate of the interpolated value
real (PREC), intent(in) :: x
!>  x-coordinates of data points in increasing order
real (PREC), intent(in), dimension(:) :: xp
!>  y-coordinates of the data points, same length as xp.
real (PREC), intent(in), dimension(:) :: fp
!>  If true, extrapolate function value at point x if x < xp(1) or
!>  x > x(size(x)).
logical, intent(in) :: ext
!>  Value to return if x < xp(1)
real (PREC), intent(in) :: left
!>  Value to return if x > xp(size(xp))
real (PREC), intent(in) :: right
!>  Interpolated value
real (PREC), intent(out) :: fx

integer (INTSIZE) :: ilb, iub
real (PREC) :: slope

if (.not. ext) then
    if (x < xp(1)) then
        fx = left
        return
    else if (x > xp(size(xp))) then
        fx = right
        return
    end if
end if

! default: either there is a bracketing interval, or we extrapolate
! values outside of domain
ilb = interp_find (x, xp)
iub = ilb + 1
slope = (x - xp(ilb)) / (xp(iub) - xp(ilb))
fx = (1.0_PREC-slope) * fp(ilb) + slope * fp(iub)
