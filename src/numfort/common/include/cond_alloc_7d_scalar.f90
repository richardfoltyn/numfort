
! Implementation for 1d-arrays and scalar SOURCE arguments

integer, intent(in), dimension(:) :: shp
    !*  Desired shape of allocated input array.
integer, intent(out), optional :: stat
    !*  Status flag. If present, it is set to -1 if no allocation is required
    !   (because ARR already is allocated and has desired shape),
    !   to -2 if an invalid input argument is encountered, and
    !   contains the value of STAT returned by the underlying ALLOCATE otherwise.
    !   In particular, if ARR was successfully allocated STAT has value 0.

integer, parameter :: NDIM = 7
integer :: lstat

lstat = COND_ALLOC_STAT_IS_ALLOCATED

if (size(shp) /= NDIM) then
    lstat = COND_ALLOC_STAT_INVALID_ARG
    goto 100
end if

if (allocated(arr)) then
    if (.not. has_shape (arr, shp)) deallocate (arr)
end if

if (.not. allocated(arr)) then
    if (present(source)) then
        allocate (arr(shp(1),shp(2),shp(3),shp(4),shp(5),shp(6),shp(7)), &
            source=source, stat=lstat)
    else
        allocate (arr(shp(1),shp(2),shp(3),shp(4),shp(5),shp(6),shp(7)), &
            stat=lstat)
    end if
end if

100 continue
    if (present(stat)) stat = lstat
