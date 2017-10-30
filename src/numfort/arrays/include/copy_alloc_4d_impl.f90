
integer, dimension(ND) :: shp

if (.not. shape_equal (src, dst) .and. allocated(dst)) then
    deallocate (dst)
end if

if (present(src)) then
    shp = shape(src)
    if (.not. allocated(dst)) then
        allocate (dst(shp(1),shp(2),shp(3),shp(4)), source=src)
    else
        dst(:,:,:,:) = src
    end if
end if
