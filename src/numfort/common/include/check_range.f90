character (*), intent(in), optional :: name
type (status_t), intent(inout) :: status
character (*), intent(out), optional :: msg

logical :: lb_fail, ub_fail

call clear_status (status, msg)


lb_fail = .false.
if (present(lb) .and. present(val)) then
    lb_fail = val < lb
end if

ub_fail = .false.
if (present(ub) .and. present(val)) then
    ub_fail = val > ub
end if

status = NF_STATUS_OK
if (lb_fail .or. ub_fail) status = NF_STATUS_INVALID_ARG

if (present(msg) .and. (lb_fail .or. ub_fail)) then
    if (present(name)) then
        msg = "Invalid argument '" // name // ": "
    else
        msg = "Invalid argument: "
    end if
    if (lb_fail) then
        msg = trim(msg) // "value violates lower bound"
    else
        msg = trim(msg) // "value violates upper bound"
    end if
end if

