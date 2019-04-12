
    ! Implementation of WRITE_FIXED

    character (*), intent(in) :: path
    character (*), intent(in), optional :: fmt
    character (*), intent(in), optional :: transform
        !*  Ignored, present only for API compatibility with other READ_FIXED
        !   routines.
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: iostat, uid
    character (100) :: lmsg
    type (status_t) :: lstatus
    integer :: n, k
    character (:), allocatable :: lfmt, fmt_field

    lstatus = NF_STATUS_IO_ERROR
    lmsg = ""

    n = size(dat, 1)

    if (present(fmt))  then
        k = len(fmt)
        if (index(fmt, '(') > 0) then
            allocate (lfmt, source=trim(fmt))
        else
            ! Caller passed element-specific format, embed this in a format
            ! string that creates the default number of columns
            allocate (character (k+20) :: lfmt)
            allocate (character (k) :: fmt_field)
            call format_strip_parenthesis (fmt, fmt_field)
            write (lfmt, '("(", i0, "(", a, "))")') n, trim(fmt_field)

            deallocate (fmt_field)
        end if
    else
        allocate (character (100) :: lfmt)
        write (lfmt, '("(*(", i0, "(g16.8),:,/))")') n
    end if

    open (newunit=uid, file=path, action='write', &
        access='sequential', iostat=iostat, iomsg=lmsg, err=100)

    write (unit=uid, fmt=lfmt, iomsg=lmsg, iostat=iostat, err=50) dat

    lstatus = NF_STATUS_OK

50  continue

    close (unit=uid)

100 continue
    if (present(msg)) lmsg = msg
    if (present(status)) status = lstatus
