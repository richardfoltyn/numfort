
logical :: res
integer, dimension(ND) :: shp1, shp2

shp1 = shape(arr1)
shp2 = shape(arr2)

res = all(shp1==shp2)
