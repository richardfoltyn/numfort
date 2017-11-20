
real (PREC), intent(in), dimension(:) :: x
    !*  1-D input array.
real (PREC), intent(out), dimension(:,:) :: xp
    !*  Matrix of powers of X, such that each row corresponds to one element
    !   of X and each column to an exponent in {0,1,...,NP-1}, where
    !   NP=SIZE(XP,2).
logical, intent(in), optional :: increasing
    !*  If present and .TRUE., sort powers in increasing order, ie.
    !   first column contains X**0.0 = 1.0.

logical :: lincr
integer :: nx, np, i, step, ifrom, ito

lincr = .false.
if (present(increasing)) lincr = increasing

nx = min(size(x), size(xp,1))
np = size(xp,2)

if (nx == 0 .or. np == 0) return

if (lincr) then
    step = 1
    ifrom = 2
    ito = np
    xp(:,1) = 1
else
    step = -1
    ifrom = np-1
    ito = 1
    xp(:,np) = 1
end if

do i = ifrom, ito, step
    xp(1:nx,i) = xp(1:nx,i-step) * x(1:nx)
end do
