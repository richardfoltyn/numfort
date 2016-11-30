integer, intent(out), dimension(:) :: rnk
    !!  Array of indices that sort arr.
character (*), intent(in), optional :: algorithm
    !!  Sorting algorithm to use. Valid values are: "mergesort"
integer, intent(out), optional :: status
    !!  Status flag.

if (size(arr) /= size(rnk)) then
    if (present(status)) status = STATUS_INVALID_INPUT
    return
end if

call mrgrnk (arr, rnk)

if (present(status)) status = STATUS_OK
