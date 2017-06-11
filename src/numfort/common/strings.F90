module numfort_common_strings
    !*  Module containing helper functions for working with character
    !   data types. We don't implement a string data type since for numfort's
    !   purposes a few basic routines are sufficient.

    use, intrinsic :: iso_fortran_env
    implicit none

    private

    public :: lower

    integer, parameter :: ASCII_LOWER_A = iachar('a')
    integer, parameter :: ASCII_LOWER_Z = iachar('z')

    integer, parameter :: ASCII_UPPER_A = iachar('A')
    integer, parameter :: ASCII_UPPER_Z = iachar('Z')


contains

subroutine lower (s)
    !*  LOWER transforms character strings to lower case in place, leaving
    !   any characters outside of the range A-Z (in ASCII encoding) unchanged.

    character (*), intent(in out) :: s
        !*  Character string that to be transformed to lower case. On exit,
        !   contains transformed string.

    integer :: i, j, offset

    offset = ASCII_LOWER_A - ASCII_UPPER_A

    do i = 1, len(s)
        j = iachar (s(i:i))
        if (j >= ASCII_UPPER_A .and. j <= ASCII_UPPER_Z) then
            s(i:i) = achar (offset + j)
        end if
    end do

end subroutine

end module
