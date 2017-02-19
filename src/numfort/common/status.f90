module numfort_common_status

    use numfort_common_enums
    implicit none
    private

    public :: decode_status
contains

pure subroutine decode_status (status, x, n)
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

    if (size(x) >= 1 .and. ior(status, 0) == 0) then
        n = 1
        x(1) = 0
    else
        n = 0
        do i = 1, size(x)
            pattern = ishft(1, i-1)
            if (iand(status, pattern) == pattern) then
                n = n + 1
                x(n) = i
            end if
        end do
    end if
end subroutine
end module
