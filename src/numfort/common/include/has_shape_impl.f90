! Implemenentation of HAS_SHAPE that contains code that is independent
! of rank, type and kind.

integer, intent(in), dimension(:) :: shp
    !*  Shape which ARR shape should be compared to.
logical :: res
    !*  True if given array has the desired shape and false otherwise.

integer, dimension(NDIM) :: arr_shp
logical, dimension(NDIM) :: dim_equal

arr_shp = shape(arr)
dim_equal = (arr_shp == shp)
res = all(dim_equal)
