


module numfort_common_print_helpers

    use, intrinsic :: iso_fortran_env

    use numfort_common_enums

    implicit none

    private


    public :: print_value
    public :: print_section
    public :: print_msg


    interface print_value
        procedure print_int, print_real64, print_1d_real64, print_2d_real64
    end interface

    contains


subroutine print_prefix (prefix, buffer)
    character (*), intent(in), optional :: prefix
    character (*), intent(inout) :: buffer

    integer :: n

    if (.not. present(prefix)) return

    n = len_trim (prefix)
    write (buffer, '("(a", i0, ")")') n + 2
    write (ERROR_UNIT, buffer, advance='NO') trim(prefix) // ": "

end subroutine



subroutine print_indent (indent, buffer)
    integer, intent(in), optional :: indent
    character (*), intent(inout) :: buffer

    if (.not. present(indent)) return

    if (indent > 0) then
        write (buffer, '("(tr", i0, ")")') indent
        write (ERROR_UNIT, buffer, advance='NO')
    end if

end subroutine



subroutine print_msg (msg, iprint, min_level, prefix, indent)
    character (*), intent(in) :: msg
    integer (NF_ENUM_KIND), intent(in) :: iprint
    integer (NF_ENUM_KIND), intent(in) :: min_level
    character (*), intent(in), optional :: prefix
    integer, intent(in), optional :: indent

    integer :: n
    character (:), allocatable :: buffer

    if (iprint < min_level) return

    allocate (character (20) :: buffer)

    call print_prefix (prefix, buffer)
    call print_indent (indent, buffer)

    ! print message
    n = len_trim (msg)
    write (buffer, '("(a", i0, ")")') n
    write (ERROR_UNIT, buffer) trim(msg)

    deallocate (buffer)

end subroutine



subroutine print_int (msg, value, iprint, min_level, prefix, indent, fmt)
    character (*), intent(in) :: msg
    integer, intent(in) :: value
    integer (NF_ENUM_KIND), intent(in) :: iprint
    integer (NF_ENUM_KIND), intent(in) :: min_level
    character (*), intent(in), optional :: prefix
    integer, intent(in), optional :: indent
    character (*), intent(in), optional :: fmt

    character (:), allocatable :: buffer, lfmt
    integer :: n

    if (iprint < min_level) return

    n = 0
    if (present(fmt)) n = len_trim (fmt)
    allocate (character (20+n) :: buffer)

    call print_prefix (prefix, buffer)
    call print_indent (indent, buffer)

    if (present(fmt)) then
        allocate (lfmt, source=trim(fmt))
    else
        allocate (lfmt, source='i0')
    end if

    n = len (msg)
    write (buffer, '("(a", i0, ", ", a10, ")")') n, lfmt
    write (ERROR_UNIT, buffer) msg, value

    deallocate (buffer, lfmt)

end subroutine



subroutine print_real64 (msg, value, iprint, min_level, prefix, indent, fmt)
    integer, parameter :: PREC = real64
    character (*), intent(in) :: msg
    real (PREC), intent(in) :: value
    integer (NF_ENUM_KIND), intent(in) :: iprint
    integer (NF_ENUM_KIND), intent(in) :: min_level
    character (*), intent(in), optional :: prefix
    integer, intent(in), optional :: indent
    character (*), intent(in), optional :: fmt

    character (:), allocatable :: buffer, lfmt
    integer :: n

    if (iprint < min_level) return

    n = 0
    if (present(fmt)) n = len_trim (fmt)
    allocate (character (30+n) :: buffer)

    call print_prefix (prefix, buffer)
    call print_indent (indent, buffer)

    if (present(fmt)) then
        allocate (lfmt, source=trim(fmt))
    else
        allocate (lfmt, source='g24.8')
    end if

    n = len (msg)

    write (buffer, '("(a", i0, ", ", a20, ")")') n, lfmt
    write (ERROR_UNIT, buffer) msg, value

    deallocate (buffer, lfmt)

end subroutine



subroutine print_1d_real64 (msg, value, iprint, min_level, prefix, indent, fmt)
    integer, parameter :: PREC = real64
    character (*), intent(in) :: msg
    real (PREC), intent(in), dimension(:) :: value
    integer (NF_ENUM_KIND), intent(in) :: iprint
    integer (NF_ENUM_KIND), intent(in) :: min_level
    character (*), intent(in), optional :: prefix
    integer, intent(in), optional :: indent
    character (*), intent(in), optional :: fmt

    character (:), allocatable :: buffer, lfmt
    integer :: n

    if (iprint < min_level) return

    n = 0
    if (present(fmt)) n = len_trim (fmt)
    allocate (character (40+n) :: buffer)

    call print_prefix (prefix, buffer)
    call print_indent (indent, buffer)

    if (present(fmt)) then
        allocate (lfmt, source=trim(fmt))
    else
        allocate (lfmt, source='g24.8')
    end if

    n = len (msg)

    write (buffer, '("(a", i0, ", *(", a20, ",:,", a4, "))")') n, lfmt, '", "'
    write (ERROR_UNIT, buffer) msg, value

    deallocate (buffer, lfmt)

end subroutine



subroutine print_2d_real64 (msg, value, iprint, min_level, prefix, indent, fmt)
    integer, parameter :: PREC = real64
    character (*), intent(in) :: msg
    real (PREC), intent(in), dimension(:,:) :: value
    integer (NF_ENUM_KIND), intent(in) :: iprint
    integer (NF_ENUM_KIND), intent(in) :: min_level
    character (*), intent(in), optional :: prefix
    integer, intent(in), optional :: indent
    character (*), intent(in), optional :: fmt

    character (:), allocatable :: buffer, lfmt, buffer2
    real (PREC), dimension(:), allocatable :: rwork
    integer :: n, i

    if (iprint < min_level) return

    n = 0
    if (present(fmt)) n = len_trim (fmt)
    allocate (character (40+n) :: buffer)
    ! Create separate buffer for prefix/indentation routines
    allocate (character (10) :: buffer2)

    if (present(fmt)) then
        allocate (lfmt, source=trim(fmt))
    else
        allocate (lfmt, source='g24.8')
    end if

    n = len (msg)

    write (buffer, '("(a", i0, ", *(", a20, ",:,", a4, "))")') n, lfmt, '", "'

    allocate (rwork(size(value, 2)))

    do i = 1, size(value, 1)
        call print_prefix (prefix, buffer2)
        call print_indent (indent, buffer2)

        rwork(:) = value(i,:)

        if (i == 1) then
            write (ERROR_UNIT, buffer) msg, rwork
        else
            write (ERROR_UNIT, buffer) "", rwork
        end if
    end do

    deallocate (buffer, buffer2, lfmt)

end subroutine



subroutine print_section (msg, iprint, min_level, prefix)
    character (*), intent(in) :: msg
    integer (NF_ENUM_KIND), intent(in) :: iprint
    integer (NF_ENUM_KIND), intent(in) :: min_level
    character (*), intent(in), optional :: prefix

    character (20) :: buffer
    integer :: n, nprefix

    if (iprint >= min_level) then

        n = len_trim (msg)

        if (present(prefix)) then
            nprefix = len_trim (prefix)
            write (buffer, '("(a", i0, ")")') nprefix + n + 2
            write (ERROR_UNIT, buffer) trim(prefix) // ": " // repeat("=", n)
            write (ERROR_UNIT, buffer) trim(prefix) // ": " // trim(msg)
            write (ERROR_UNIT, buffer) trim(prefix) // ": " // repeat("=", n)
        else
            write (buffer, '("(a", i0, ")")') n
            write (ERROR_UNIT, buffer) repeat("=", n)
            write (ERROR_UNIT, buffer) trim(msg)
            write (ERROR_UNIT, buffer) repeat("=", n)
        end if

    end if

end subroutine



end module
