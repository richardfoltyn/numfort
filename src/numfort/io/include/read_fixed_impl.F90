
    ! Implementation of READ_FIXED

    character (*), intent(in) :: path
    character (*), intent(in) :: fmt
    character (*), intent(in), optional :: transform
        !*  Ignored, present only for API compatibility with other READ_FIXED
        !   routines.
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: iostat, uid
    character (100) :: iomsg
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    if (present(msg)) msg = ""

    open (newunit=uid, file=path, status='old', action='read', &
        access='sequential', iostat=iostat, iomsg=iomsg)

    if (iostat /= 0) then
        lstatus = NF_STATUS_IO_ERROR
        if (present(msg)) msg = iomsg
        goto 100
    end if

    read (unit=uid, fmt=fmt, iomsg=iomsg, iostat=iostat) dat

    if (iostat /= 0) then
        lstatus = NF_STATUS_IO_ERROR
        if (present(msg)) msg = iomsg
        goto 50
    end if

50  continue

    close (unit=uid)

100 continue

    if (present(status)) status = lstatus
