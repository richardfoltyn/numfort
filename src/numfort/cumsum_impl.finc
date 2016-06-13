integer, intent(in), dimension(:) :: shp
integer, intent(in) :: axis

integer, dimension(size(shp)) :: strides
integer :: i, k, rnk, offset
integer :: nstrides_h, stride_h, stride_axis

integer :: ilb, iub, ilb_prev, iub_prev

rnk = size(shp)

strides = 1
do i = 2, rnk
    strides(i) = strides(i-1) * shp(i-1)
end do

nstrides_h = product(shp(axis + 1:))

stride_h = 0

if (axis < rnk) then
    stride_h = strides(axis + 1)
end if

res = 0
stride_axis = strides(axis)

if (axis == 1) then
    do k = 1, nstrides_h
        offset = (k-1) * stride_h + 1
        res(offset) = x(offset)
        do i = offset + 1, offset + shp(axis) - 1
            res(i) = res(i - 1) + x(i)
        end do
    end do
else
    ! iterate over dimenions that come AFTER axis
    do k = 1, nstrides_h
        offset = (k-1) * stride_h + 1
        ! copy over first block at index=1 on for relevant axis
        ilb = offset
        iub = offset + stride_axis - 1
        res(ilb:iub) = x(ilb:iub)
        ! iterate over axis which we apply reduction operation over
        do i = 2, shp(axis)
            ilb = offset + (i-1) * stride_axis
            iub = ilb + stride_axis - 1
            ilb_prev = ilb - stride_axis
            iub_prev = ilb - 1
            res(ilb:iub) = res(ilb_prev:iub_prev) + x(ilb:iub)
        end do
    end do
end if
