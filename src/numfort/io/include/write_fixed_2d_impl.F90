
    ! Implementation of WRITE_FIXED

    character (*), intent(in) :: path
    character (*), intent(in) :: fmt
    character (*), intent(in), optional :: transform
        !   routines.
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: iostat, n, uid, i
    integer (NF_ENUM_KIND) :: itransform
    character (100) :: iomsg
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    if (present(msg)) msg = ""

    itransform = parse_transform (transform)
    if (itransform == 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    open (newunit=uid, file=path, action='write', &
        access='sequential', iostat=iostat, iomsg=iomsg)

    if (iostat /= 0) then
        lstatus = NF_STATUS_IO_ERROR
        if (present(msg)) msg = iomsg
        goto 100
    end if

    select case (itransform)
    case (NF_IO_TRANSFORM_NONE)
        ! Write data such that one row (line) in the data file corresponds
        ! to one row in the output array in the same order
        n = size(dat, 1)
        do i = 1, n
            write (unit=uid, fmt=fmt, iomsg=iomsg, iostat=iostat) dat(i,:)
            if (iostat /= 0) then
                lstatus = NF_STATUS_IO_ERROR
                if (present(msg)) msg = iomsg
                goto 50
            end if
        end do

    case (NF_IO_TRANSFORM_TRANSPOSE)
        ! Write data such that one row (line) in the data file corresponds
        ! to one column in the output array
        n = size(dat, 1)
        do i = 1, n
            write (unit=uid, fmt=fmt, iomsg=iomsg, iostat=iostat) dat(:,i)
            if (iostat /= 0) then
                lstatus = NF_STATUS_IO_ERROR
                if (present(msg)) msg = iomsg
                goto 50
            end if
        end do

    case (NF_IO_TRANSFORM_FLATTEN)
        ! Read in one go in column-major order
        write (unit=uid, fmt=fmt, iomsg=iomsg, iostat=iostat) dat
        if (iostat /= 0) then
            lstatus = NF_STATUS_IO_ERROR
            if (present(msg)) msg = iomsg
            goto 50
        end if

    end select

50  continue

    close (unit=uid)

100 continue

    if (present(status)) status = lstatus
