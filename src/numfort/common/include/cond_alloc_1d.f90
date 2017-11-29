
! Implementation for 1d-arrays and array SOURCE arguments

integer, intent(out), optional :: stat
    !*  Status flag. If present, it is set to -1 if no allocation is required
    !   (because ARR already is allocated and has desired shape), and
    !   contains the value of STAT returned by the underlying ALLOCATE otherwise.
    !   In particular, if ARR was successfully allocated STAT has value 0.

integer, dimension(1) :: shp
integer :: lstat

lstat = -1
shp = shape(source)

if (allocated(arr)) then
    if (.not. shape_equal (arr, source)) deallocate (arr)
end if

if (.not. allocated(arr)) then
    allocate (arr(shp(1)), source=source, stat=lstat)
end if

if (present(stat)) stat = lstat
