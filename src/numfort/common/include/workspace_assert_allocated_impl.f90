! Type-independent part of workspace_assert_allocated

intent(in out) :: arr
integer, intent(in out) :: n_attr
integer, intent(in), optional :: n

allocatable :: arr, tmp

if (present(n)) then
    if (n > 0) then
        if (.not. allocated(arr)) then
            allocate (arr(n))
            n_attr = n
        else if (size(arr) < n) then
            allocate (tmp(n))
            tmp(1:size(arr)) =  arr
            call move_alloc (tmp, arr)
            n_attr = n
        end if
    end if
end if
