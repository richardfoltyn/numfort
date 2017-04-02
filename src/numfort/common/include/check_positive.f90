character (*), intent(in), optional :: name
type (status_t), intent(in out) :: status
character (*), intent(out), optional :: msg

call clear_status (status, msg)

if (present(val)) then
    if (val <= 0) then
        call status_set (status, NF_STATUS_INVALID_ARG)
        if (present(msg)) then
            if (present(name)) then
                msg = "Invalid argument '" // name // "': positive real number required"
            else
                msg = "Invalid argument: positive real number required"
            end if
        end if
    end if
end if
