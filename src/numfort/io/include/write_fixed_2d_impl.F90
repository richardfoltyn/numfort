
    ! Implementation of WRITE_FIXED

    character (*), intent(in) :: path
    character (*), intent(in) :: fmt
    character (*), intent(in), optional :: transform
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: iostat, n, uid, i, nlines, ncol
    integer (NF_ENUM_KIND) :: itransform
    character (100) :: lmsg
    character (:), allocatable :: lfmt, fmt_field
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    lmsg = ""

    itransform = parse_transform (transform)
    if (itransform == 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Default status for remainder of routine
    lstatus = NF_STATUS_IO_ERROR

    open (newunit=uid, file=path, action='write', &
        access='sequential', iostat=iostat, iomsg=lmsg, err=100)

    ! Format string used for TRANSFORM_NONE and TRANSFORM_TRANSPOSE
    n = len_trim (fmt)
    allocate ( character(n+20) :: lfmt)
    allocate ( character(n) :: fmt_field)
    ! Strip original format of all parenthesis to obtain the format for
    ! individual values within a record without parenthesis.
    call format_strip_parenthesis (fmt, fmt_field)

    select case (itransform)
    case (NF_IO_TRANSFORM_NONE)
        ! Write data such that one row (line) in the data file corresponds
        ! to one column in the input array, thus preserving the linear order
        ! as found in the data file also in memory

        ! Number of columns IN THE DATA FILE
        ncol = size(dat, 1)
        write (lfmt, '("(*(", i0, "(", a, "),:,/))")') ncol, trim(fmt_field)

        write (unit=uid, fmt=lfmt, iomsg=lmsg, iostat=iostat, err=50) dat

    case (NF_IO_TRANSFORM_TRANSPOSE)
        ! Write data such that one row (line) in the data file corresponds to
        ! one row in the output array. This effectively transposes the data,
        ! and the linear order is no longer preserved as Fortran uses
        ! column-major ordering.

        ! Number of lines and columns IN THE DATA FILE
        nlines = size(dat, 1)
        ncol = size(dat, 2)
        write (lfmt, '("(", i0, "(", a, "))")') ncol, trim(fmt_field)

        ! Use contiguous buffer for I/O calls to prevent temp. array creation
        allocate (buf(ncol))

        do i = 1, nlines
            buf(:) = dat(i,:)
            write (unit=uid, fmt=lfmt, iomsg=lmsg, iostat=iostat, err=50)  buf
        end do

        deallocate (buf)

    case (NF_IO_TRANSFORM_FORMAT)
        ! Use unmodified user-provided format specifier, write stuff in
        ! column-major order.
        write (unit=uid, fmt=fmt, iomsg=lmsg, iostat=iostat, err=50) dat

    end select

    lstatus = NF_STATUS_OK

50  continue

    close (unit=uid)

100 continue
    if (present(msg)) msg = lmsg
    if (present(status)) status = lstatus
