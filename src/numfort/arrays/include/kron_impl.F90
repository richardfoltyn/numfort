

integer :: i1, j1, j2, ic1, ic, ir1
integer :: m1, m2, n1, n2, m, n
integer :: shp(2)

m1 = size(mat1, 1)
n1 = size(mat1, 2)
m2 = size(mat2, 1)
n2 = size(mat2, 2)

m = m1 * m2
n = n1 * n2

shp(1) = m
shp(2) = n
if (.not. has_shape(res, shp)) then
    return
end if


do j1 = 1, n1
    ! column offset on output matrix due to matrix 1
    ic1 = (j1-1) * n2
    do j2 = 1, n2
        ! column index on output matrix
        ic = ic1 + j2

        do i1 = 1, m1
            ! Row offset on output matrix
            ir1 = (i1-1) * m2
            res(ir1+1:ir1+m2,ic) = mat1(i1,j1) * mat2(1:m2,j2)
        end do
    end do
end do

