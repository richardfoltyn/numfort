! Implementation for 1D linear interpolation
    
real (PREC), intent(in), dimension(:) :: x, xp, fp
logical, intent(in), optional :: ext
real (PREC), intent(in), optional :: left, right
real (PREC), intent(out), dimension(:) :: fx
    
integer (INTSIZE) :: n, i, ilb, iub
logical :: lext
real (PREC) :: lright, lleft, slope, xi, xp1, xpn
    
n = size (xp)
lext = .false.
lleft = fp(1)
lright = fp(n)
    
if (present(ext)) lext = ext
if (present(left)) lleft = left
if (present(right)) lright = right
    
xp1 = xp(1)
xpn = xp(n)
    
do i = 1, n
    xi = x(i)        
    if (.not. lext) then
        if (xi < xp1) then
            fx(i) = lleft
            cycle
        else if (xi > xpn) then
            fx(i) = lright
            cycle
        end if
    end if
            
    ! default: either there is a bracketing interval, or we extrapolate
    ! values outside of domain
    ilb = interp_find(xi, xp)
    iub = ilb + 1
    slope = (xi - xp(ilb)) / (xp(iub) - xp(ilb))
    fx(i) = (1-slope) * fp(ilb) + slope * fp(iub)
        
end do