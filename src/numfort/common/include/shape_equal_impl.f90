
logical :: res
logical, dimension(NDIM) :: all_equal
integer, dimension(NDIM) :: shp1, shp2

res = .false.
if (.not. present(arr2)) return

shp1 = shape(arr1)
shp2 = shape(arr2)

all_equal = (shp1 == shp2)
res = all(all_equal)
