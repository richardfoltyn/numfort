! Implementation of arange() subroutine

real (PREC), intent(out), dimension(:) :: x
real (PREC), intent(in), optional :: ifrom
real (PREC), intent(in), optional :: step
    
real (PREC) :: val, lstep, lifrom
integer :: i

lstep = 1.0_PREC
lifrom = 1.0_PREC
if (present(step)) lstep = step
if (present(ifrom)) lifrom = ifrom
    
val = lifrom
do i = 1, size(x)
    x(i) = val
    val = val + lstep
end do

