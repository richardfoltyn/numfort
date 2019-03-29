
    ! Implementation of WRITE_FIXED

    character (*), intent(in) :: path
    character (*), intent(in) :: fmt
    character (*), intent(in), optional :: transform
        !*  Ignored, present only for API compatibility with other READ_FIXED
        !   routines.
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: iostat, uid
    character (100) :: lmsg
    type (status_t) :: lstatus

    lstatus = NF_STATUS_IO_ERROR
    lmsg = ""

    open (newunit=uid, file=path, action='write', &
        access='sequential', iostat=iostat, iomsg=lmsg, err=100)

    write (unit=uid, fmt=fmt, iomsg=lmsg, iostat=iostat, err=50) dat

    lstatus = NF_STATUS_OK

50  continue

    close (unit=uid)

100 continue
    if (present(msg)) lmsg = msg
    if (present(status)) status = lstatus
