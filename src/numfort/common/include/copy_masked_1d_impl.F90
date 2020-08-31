
real (PREC), intent(in), dimension(:), contiguous :: src
real (PREC), intent(out), dimension(:), contiguous :: dst
logical, intent(in), dimension(:), contiguous :: mask
type (status_t), intent(out), optional :: status

type (status_t) :: lstatus
integer :: n

lstatus = NF_STATUS_OK

if (size(src) /= size(mask)) then
    lstatus = NF_STATUS_INVALID_ARG
    goto 100
end if

n = count (mask)

if (size(dst) < n) then
    lstatus = NF_STATUS_INVALID_ARG
    goto 100
end if

dst(1:n) = pack (src, mask)

100 continue

if (present(status)) status = lstatus
