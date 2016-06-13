integer, optional :: axis

integer :: laxis
integer, dimension(size(shape(x))) :: shp

laxis = 1
if (present(axis)) laxis = axis

shp = shape(x)
ptr_x(1:size(x)) => x
ptr_res(1:size(res)) => res

call cumsum_nd_real64(ptr_x, shp, ptr_res, laxis)
