integer, intent(out), dimension(:) :: rnk
    !!  Array of indices that sort arr.
character (*), intent(in), optional :: algorithm
    !!  Sorting algorithm to use. Valid values are: "mergesort"
type (status_t), intent(out), optional :: status
    !!  Status flag.

if (size(arr) /= size(rnk)) then
    if (present(status)) then
        call status_set (status, NF_STATUS_INVALID_ARG)
    end if
    return
end if

call mrgrnk (arr, rnk)

if (present(status)) then
    call status_set (status, NF_STATUS_OK)
end if
