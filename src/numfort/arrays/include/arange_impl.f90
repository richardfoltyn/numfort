! Implementation of arange() subroutine

integer (INTSIZE), intent(out), dimension(:) :: x
integer (INTSIZE), intent(in), optional :: ifrom
integer (INTSIZE), intent(in), optional :: step
    
integer (INTSIZE) :: i, val, lstep, lifrom

lstep = 1
lifrom = 1
if (present(step)) lstep = step
if (present(ifrom)) lifrom = ifrom
    
val = lifrom
do i = 1, size(x)
    x(i) = val
    val = val + lstep
end do
    