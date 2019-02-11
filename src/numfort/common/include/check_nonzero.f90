character (*), intent(in), optional :: name
type (status_t), intent(inout) :: status
character (*), intent(out), optional :: msg

call clear_status (status, msg)

if (present(val)) then
    if (val == 0) then
        status = NF_STATUS_INVALID_ARG
        if (present(msg)) then
            if (present(name)) then
                msg = "Invalid argument '" // name // "': non-zero value required"
            else
                msg = "Invalid argument: non-zero value required"
            end if
        end if
    end if
end if
