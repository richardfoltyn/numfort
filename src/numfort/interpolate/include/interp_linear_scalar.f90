real (PREC), intent(in) :: x
real (PREC), intent(in), dimension(:) :: xp, fp
logical, intent(in), optional :: ext
real (PREC), intent(in), optional :: left, right
real (PREC), intent(out) :: fx

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
