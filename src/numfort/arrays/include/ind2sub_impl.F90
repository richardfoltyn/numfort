integer (INTSIZE), intent(in), dimension(:) :: shp
    !!  Shape of array to use to unravel flat indices.
integer (INTSIZE), intent(in), dimension(:) :: lin_indices
    !!  Array of flat indices to convert.
integer (INTSIZE), intent(out), dimension(:,:) :: sub_indices
    !!  Array to store unraveled coordinate tuples.

integer (INTSIZE), dimension(size(shp)) :: stride
integer (INTSIZE) :: i, j, n, d, rnk, stride_j, rem_i, s_ij
integer (INTSIZE), dimension(:), allocatable :: rem

rnk = size(shp)

stride = 1
do i = 2, rnk
    stride(i) = stride(i-1) * shp(i-1)
end do

n = size(sub_indices, 1)
d = size(sub_indices, 2)

! automatically deallocated on procedure exit
allocate (rem(n), source=lin_indices)

rem = rem - 1

do j = d, 2, -1
    stride_j = stride(j)
    do i = 1, n
        ! current remainder for i
        rem_i = rem(i)
        ! keep only integer part of expression; no need to apply floor()
        s_ij = rem_i / stride_j
        ! adjust remainder for processed dimension
        rem_i = rem_i - s_ij * stride_j
        ! adjust index to start at 1
        sub_indices(i, j) = s_ij + 1
        rem(i) = rem_i
    end do
end do

! index on first dimension is just the remainder after all higher dimensions
! have been subtracted out
sub_indices(:, 1) = rem + 1
