integer (INTSIZE), intent(out), dimension(:) :: rnk
    !*  Array of indices that sort arr.
character (*), intent(in), optional :: algorithm
    !*  Sorting algorithm to use. Valid values are: "mergesort"
integer (NF_ENUM_KIND), intent(out), optional :: status
    !*  Status flag.

if (size(arr) /= size(rnk)) then
    if (present(status)) status = NF_STATUS_INVALID_ARG
    return
end if

call mrgrnk (arr, rnk)

if (present(status)) status = NF_STATUS_OK
