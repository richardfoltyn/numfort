

module numfort_io_common

    use, intrinsic :: iso_fortran_env

    use numfort_common

    implicit none
    private


    public :: data_chunk_real32
    public :: data_chunk_real64
    public :: get_temp_directory
    public :: format_get_field_width
    public :: format_strip_parenthesis

    integer (NF_ENUM_KIND), public, parameter :: NF_IO_TRANSFORM_NONE = 1
    integer (NF_ENUM_KIND), public, parameter :: NF_IO_TRANSFORM_TRANSPOSE = 2
    integer (NF_ENUM_KIND), public, parameter :: NF_IO_TRANSFORM_FORMAT = 4


    type :: data_chunk_real32
        type (data_chunk_real32), pointer :: ptr_next => NULL ()
        real (real32), dimension(:,:), allocatable :: dat
    end type

    type :: data_chunk_real64
        type (data_chunk_real64), pointer :: ptr_next => NULL ()
        real (real64), dimension(:,:), allocatable :: dat
    end type


    interface format_get_field_width
        procedure format_get_field_width_real32, format_get_field_width_real64
    end interface

    contains



function format_get_field_width_real32 (fmt, dummy) result(res)
    !*  FORMAT_GET_FIELD_WIDTH attempts to determine the field implied
    !   by a format specification.
    integer, parameter :: PREC = real32
    character (*), intent(in) :: fmt
    real (PREC), intent(in), value :: dummy
        !*  Dummy argument required to correctly determine the specific
        !   routine for a generic interface.
    integer :: res

    integer, parameter :: BUF_LEN = 1024
    character (BUF_LEN) :: buf
    integer :: iostat
    character (100) :: iomsg

    write (buf, fmt, iostat=iostat, iomsg=iomsg) 1.0_PREC/9.0_PREC

    if (iostat /= 0) then
        res = -1
        return
    end if

    res = len_trim (buf)

end function



function format_get_field_width_real64 (fmt, dummy) result(res)
    !*  FORMAT_GET_FIELD_WIDTH attempts to determine the field implied
    !   by a format specification.
    integer, parameter :: PREC = real64
    character (*), intent(in) :: fmt
    real (PREC), intent(in), value :: dummy
        !*  Dummy argument required to correctly determine the specific
        !   routine for a generic interface.
    integer :: res

    integer, parameter :: BUF_LEN = 1024
    character (BUF_LEN) :: buf
    integer :: iostat
    character (100) :: iomsg

    write (buf, fmt, iostat=iostat, iomsg=iomsg) 1.0_PREC/9.0_PREC

    if (iostat /= 0) then
        res = -1
        return
    end if

    res = len_trim (buf)

end function



subroutine format_strip_parenthesis (fmt, fmt_out)
    !*  FORMAT_STRIP_PARENTHESIS returns the innermost format string
    !   contains in (a set of possibly nested) parenthesis.
    !   For example, if FMT = '(*(fW.D))', the routine will return
    !   FMT_OUT = 'fW.D'
    character (*), intent(in) :: fmt
    character (*), intent(out) :: fmt_out

    integer :: i

    fmt_out = fmt

    do while (.true.)
        i = index(fmt_out, '(')
        if (i <= 0) exit
        fmt_out = fmt_out(i+1:len(fmt_out))
    end do

    i = index (fmt_out, ')')
    if (i > 0) then
        fmt_out = fmt_out(1:i-1)
    end if

end subroutine



subroutine get_temp_directory (path, status)
    !*  GET_TEMP_DIRECTORY attempts to return the path to a system-specific
    !   temporary directory. No attempt is made to ensure that this
    !   directory actually exists or is writable by the calling process,
    !   so this should only be used internally (e.g. for unit tests).
    character (*), intent(inout) :: path
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: length

    lstatus = NF_STATUS_OK

    ! Look for TMP first
    call get_environment_variable ("TMP", path, length=length)
    if (length > 0) then
        if (len(path) < length) then
            lstatus = NF_STATUS_STORAGE_ERROR
        end if
        goto 100
    end if

    ! TMP does not exist, try TEMP
    call get_environment_variable ("TEMP", path, length=length)
    if (length > 0) then
        if (len(path) < length) then
            lstatus = NF_STATUS_STORAGE_ERROR
        end if
        goto 100
    end if

    ! TMP and TEMP do not exist, try TMPDIR
    call get_environment_variable ("TMPDIR", path, length=length)
    if (length > 0) then
        if (len(path) < length) then
            lstatus = NF_STATUS_STORAGE_ERROR
        end if
        goto 100
    end if

#if defined(_LINUX)
    path = '/tmp'
    if (len(path) < len('/tmp')) then
        lstatus = NF_STATUS_STORAGE_ERROR
    end if
#elif defined(_WIN32)
    path = 'C:/Temp'
    if (len(path) < 'C:/Temp') then
        lstatus = NF_STATUS_STORAGE_ERROR
    end if
#endif

100 continue
    if (present(status)) status = lstatus

end subroutine


end module
