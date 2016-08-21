! Implementation of arange() subroutine

integer (INTSIZE), intent(out), dimension(:) :: x
integer (INTSIZE), intent(in) :: ifrom
integer (INTSIZE), intent(in), optional :: step
    
integer (INTSIZE) :: i, val, lstep

lstep = 1
if (present(step)) lstep = step
    
val = ifrom
do i = 1, size(x)
    x(i) = val
    val = val + lstep
end do
    