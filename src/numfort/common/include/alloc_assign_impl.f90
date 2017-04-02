
integer :: n

if (present(src)) then
    n = size(src)
    ! allocate array to store
    if (allocated(dst)) then
        if (size(dst) /= n) deallocate (dst)
    end if

    if (.not. allocated(dst)) then
        allocate (dst(n), source=src)
    else
        dst(1:n) = src
    end if
end if
