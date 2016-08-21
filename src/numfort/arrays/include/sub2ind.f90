integer (INTSIZE), intent(in), dimension(:) :: shp
integer (INTSIZE), intent(in), dimension(:,:) :: sub_indices
integer (INTSIZE), intent(out), dimension(:) :: lin_indices

integer (INTSIZE), dimension(size(shape(sub_indices))) :: stride
integer (INTSIZE) :: i, j, rnk, stride_j

rnk = size(shp)

stride = 1
do i = 2, rnk
    stride(i) = stride(i-1) * shp(i-1)
end do

! the stride for the first dimension is 1, so we can just copy it over
lin_indices = sub_indices(:, 1)

! alternatively could be implemented as (M - 1_{m,d}) * s + M(:, 1),
! where M is the matrix of sub indices and s is the vector of strides.
do j = rnk, 2, -1
    stride_j = stride(j)
    do i = 1, size(lin_indices, 1)
        lin_indices(i) = lin_indices(i) + (sub_indices(i, j) - 1) * stride_j
    end do
end do
