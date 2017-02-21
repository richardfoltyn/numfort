module numfort_common_status

    use numfort_common_enums
    implicit none
    private

    public :: status_decode
    public :: status_success
contains

pure function status_success (status) result(res)
    integer (NF_ENUM_KIND), intent(in) :: status
    logical :: res
    res = (iand(status, NF_STATUS_OK) == NF_STATUS_OK)
end function

pure subroutine status_decode (status, x, n)
    !*  DECODE_STATUS disaggregates a composize status code into its
    !   components and turns their base-2 exponents.
    integer (NF_ENUM_KIND), intent(in) :: status
        !!  (Composite) status code.
    integer, dimension(:), intent(out) :: x
        !!  Array to store individual status codes. The lowest-exponent codes
        !!  up to the maximum given by size(array) are returned.
    integer, intent(out) :: n
        !!  Number of individual status codes present in composize status.
        !!  (n <= size(array)).

    integer :: i
    integer (NF_ENUM_KIND) :: pattern

    ! Since status codes are integers >= 0, initialize to something invalid.
    x = -1

    if ((size(x) >= 1) .and. (status /= NF_STATUS_UNDEFINED)) then
        n = 0
        do i = 1, min(size(x), bit_size(status))
            pattern = ishft(1, i-1)
            if (iand(status, pattern) == pattern) then
                n = n + 1
                x(n) = i-1
            end if
        end do
    end if
end subroutine
end module
