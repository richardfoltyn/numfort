! Implementation for 1D linear interpolation

real (PREC), intent(in), dimension(:) :: x, xp, fp
logical, intent(in), optional :: ext
real (PREC), intent(in), optional :: left, right
real (PREC), intent(out), dimension(:) :: fx

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
