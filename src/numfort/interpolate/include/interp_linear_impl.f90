! Implementation of linear interpolation for scalar input/return value
! and no optional arguments.
! Should be called from wrapper routines doing the optional argument handling,
! etc.

real (PREC), intent(in) :: x
real (PREC), intent(in), dimension(:) :: xp, fp
logical, intent(in) :: ext
real (PREC), intent(in) :: left, right
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
