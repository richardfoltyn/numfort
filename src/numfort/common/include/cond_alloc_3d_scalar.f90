
! Implementation for 1d-arrays and scalar SOURCE arguments

integer, intent(in), dimension(:) :: shp
    !*  Desired shape of allocated input array.
integer, intent(out), optional :: stat
    !*  Status flag. If present, it is set to -1 if no allocation is required
    !   (because ARR already is allocated and has desired shape), and
    !   contains the value of STAT returned by the underlying ALLOCATE otherwise.
    !   In particular, if ARR was successfully allocated STAT has value 0.

integer :: lstat

lstat = -1

if (allocated(arr)) then
    if (.not. has_shape (arr, shp)) deallocate (arr)
end if

if (.not. allocated(arr)) then
    if (present(source)) then
        allocate (arr(shp(1),shp(2),shp(3)), source=source, stat=lstat)
    else
        allocate (arr(shp(1),shp(2),shp(3)), stat=lstat)
    end if
end if

if (present(stat)) stat = lstat
