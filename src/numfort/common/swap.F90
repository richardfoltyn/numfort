

module numfort_common_swap
    !*  Module implements routines to swap two pointers to scalars or arrays
    !   without having to use a third temporary pointer.

    use, intrinsic :: iso_fortran_env

    implicit none

    private

    public :: swap

    interface swap
        procedure swap_1d_real32, swap_1d_real64, swap_1d_int32, swap_1d_int64
    end interface

    interface swap
        procedure swap_2d_real32, swap_2d_real64, swap_2d_int32, swap_2d_int64
    end interface

    interface swap
        procedure swap_3d_real32, swap_3d_real64, swap_3d_int32, swap_3d_int64
    end interface

    interface swap
        procedure swap_4d_real32, swap_4d_real64, swap_4d_int32, swap_4d_int64
    end interface

    interface swap
        procedure swap_real32, swap_real64
    end interface

    contains


! ------------------------------------------------------------------------------
! 1d-arrays

subroutine swap_1d_real32 (ptr1, ptr2)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout), dimension(:), pointer :: ptr1
    real (PREC), intent(inout), dimension(:), pointer :: ptr2

    real (PREC), dimension(:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_1d_real64 (ptr1, ptr2)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout), dimension(:), pointer :: ptr1
    real (PREC), intent(inout), dimension(:), pointer :: ptr2

    real (PREC), dimension(:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_1d_int32 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(inout), dimension(:), pointer :: ptr1
    integer (INTSIZE), intent(inout), dimension(:), pointer :: ptr2

    integer (INTSIZE), dimension(:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_1d_int64 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(inout), dimension(:), pointer :: ptr1
    integer (INTSIZE), intent(inout), dimension(:), pointer :: ptr2

    integer (INTSIZE), dimension(:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

! ------------------------------------------------------------------------------
! 2d-arrays

subroutine swap_2d_real32 (ptr1, ptr2)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout), dimension(:,:), pointer :: ptr1
    real (PREC), intent(inout), dimension(:,:), pointer :: ptr2

    real (PREC), dimension(:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_2d_real64 (ptr1, ptr2)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout), dimension(:,:), pointer :: ptr1
    real (PREC), intent(inout), dimension(:,:), pointer :: ptr2

    real (PREC), dimension(:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_2d_int32 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(inout), dimension(:,:), pointer :: ptr1
    integer (INTSIZE), intent(inout), dimension(:,:), pointer :: ptr2

    integer (INTSIZE), dimension(:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_2d_int64 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(inout), dimension(:,:), pointer :: ptr1
    integer (INTSIZE), intent(inout), dimension(:,:), pointer :: ptr2

    integer (INTSIZE), dimension(:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine


! ------------------------------------------------------------------------------
! 3d-arrays

subroutine swap_3d_real32 (ptr1, ptr2)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout), dimension(:,:,:), pointer :: ptr1
    real (PREC), intent(inout), dimension(:,:,:), pointer :: ptr2

    real (PREC), dimension(:,:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_3d_real64 (ptr1, ptr2)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout), dimension(:,:,:), pointer :: ptr1
    real (PREC), intent(inout), dimension(:,:,:), pointer :: ptr2

    real (PREC), dimension(:,:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_3d_int32 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(inout), dimension(:,:,:), pointer :: ptr1
    integer (INTSIZE), intent(inout), dimension(:,:,:), pointer :: ptr2

    integer (INTSIZE), dimension(:,:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_3d_int64 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(inout), dimension(:,:,:), pointer :: ptr1
    integer (INTSIZE), intent(inout), dimension(:,:,:), pointer :: ptr2

    integer (INTSIZE), dimension(:,:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine


! ------------------------------------------------------------------------------
! 4d-arrays

subroutine swap_4d_real32 (ptr1, ptr2)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout), dimension(:,:,:,:), pointer :: ptr1
    real (PREC), intent(inout), dimension(:,:,:,:), pointer :: ptr2

    real (PREC), dimension(:,:,:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_4d_real64 (ptr1, ptr2)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout), dimension(:,:,:,:), pointer :: ptr1
    real (PREC), intent(inout), dimension(:,:,:,:), pointer :: ptr2

    real (PREC), dimension(:,:,:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_4d_int32 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(inout), dimension(:,:,:,:), pointer :: ptr1
    integer (INTSIZE), intent(inout), dimension(:,:,:,:), pointer :: ptr2

    integer (INTSIZE), dimension(:,:,:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

subroutine swap_4d_int64 (ptr1, ptr2)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(inout), dimension(:,:,:,:), pointer :: ptr1
    integer (INTSIZE), intent(inout), dimension(:,:,:,:), pointer :: ptr2

    integer (INTSIZE), dimension(:,:,:,:), pointer :: ptr_tmp

    ptr_tmp => ptr1
    ptr1 => ptr2
    ptr2 => ptr_tmp
end subroutine

! ------------------------------------------------------------------------------
! Scalar routines

subroutine swap_real32 (x1, x2)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout) :: x1, x2
    real (PREC) :: tmp

    tmp = x1
    x1 = x2
    x2 = tmp
end subroutine

subroutine swap_real64 (x1, x2)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout) :: x1, x2
    real (PREC) :: tmp

    tmp = x1
    x1 = x2
    x2 = tmp
end subroutine

end module
