

module numfort_common_swap

    use, intrinsic :: iso_fortran_env

    implicit none

    private

    public :: swap

    interface swap
        procedure swap_1d_real32, swap_1d_real64, swap_1d_int32, swap_1d_int64
    end interface

    contains

subroutine swap_1d_real32 (ptr1, ptr2)
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:), pointer :: ptr1
    real (PREC), intent(in out), dimension(:), pointer :: ptr2

    real (PREC), dimension(:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_1d_real64 (ptr1, ptr2)
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:), pointer :: ptr1
    real (PREC), intent(in out), dimension(:), pointer :: ptr2

    real (PREC), dimension(:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_1d_int32 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:), pointer :: ptr1
    integer (INTSIZE), intent(in out), dimension(:), pointer :: ptr2

    integer (INTSIZE), dimension(:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_1d_int64 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:), pointer :: ptr1
    integer (INTSIZE), intent(in out), dimension(:), pointer :: ptr2

    integer (INTSIZE), dimension(:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

end module
