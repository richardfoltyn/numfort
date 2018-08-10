
! Implementation for 1d-arrays and array SOURCE arguments

integer, intent(out), optional :: stat
    !*  Status flag. If present, it is set to -1 if no allocation is required
    !   (because ARR already is allocated and has desired shape), and
    !   contains the value of STAT returned by the underlying ALLOCATE otherwise.
    !   In particular, if ARR was successfully allocated STAT has value 0.

integer, dimension(7) :: shp
integer :: lstat

lstat = -1
shp = shape(source)

if (allocated(arr)) then
    if (.not. has_shape (arr, shp)) deallocate (arr)
end if

if (.not. allocated(arr)) then
    allocate (arr(shp(1),shp(2),shp(3),shp(4),shp(5),shp(6),shp(7)), &
        source=source, stat=lstat)
end if

if (present(stat)) stat = lstat
