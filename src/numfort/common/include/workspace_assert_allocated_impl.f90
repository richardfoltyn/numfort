! Type-independent part of workspace_assert_allocated

intent(in out) :: arr
integer, intent(in), optional :: n

allocatable :: arr, tmp

if (present(n)) then
    if (.not. allocated(arr)) then
        allocate (arr(n))
    else if (size(arr) < n) then
        allocate (tmp(n))
        call move_alloc (tmp, arr)
    end if
end if
