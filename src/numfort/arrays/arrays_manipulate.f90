

module numfort_arrays_manipulate
    !*  Module contains various array manipulation routines that do not
    !   fit any more specific category.

    use, intrinsic :: iso_fortran_env

    implicit none
    private

    public :: insert


    interface insert
        procedure insert_scalar_int32, insert_scalar_real32, insert_scalar_real64
    end interface

    interface insert
        procedure insert_1d_int32, insert_1d_real32, insert_1d_real64
    end interface

    interface insert
        procedure insert_inplace_scalar_int32, insert_inplace_scalar_real32, &
            insert_inplace_scalar_real64
    end interface

    interface insert
        procedure insert_inplace_1d_int32, insert_inplace_1d_real32, &
            insert_inplace_1d_real64
    end interface

    contains


!-------------------------------------------------------------------------------
! INSERT routines

subroutine insert_1d_int32 (arr, idx, val, out)
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:) :: arr
    integer, intent(in) :: idx
    integer (INTSIZE), intent(in), dimension(:) :: val
    integer (INTSIZE), intent(out), dimension(:) :: out

    include "include/insert_impl.f90"
end subroutine

subroutine insert_inplace_1d_int32 (arr, idx, val)
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in out), dimension(:) :: arr
    integer, intent(in) :: idx
    integer (INTSIZE), intent(in), dimension(:) :: val

    include "include/insert_inplace_impl.f90"
end subroutine

subroutine insert_scalar_int32 (arr, idx, val, out)
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:) :: arr
    integer, intent(in) :: idx
    integer (INTSIZE), intent(in)  :: val
    integer (INTSIZE), intent(out), dimension(:) :: out

    integer (INTSIZE), dimension(1) :: val1d

    val1d(1) = val
    call insert (arr, idx, val1d, out)
end subroutine

subroutine insert_inplace_scalar_int32 (arr, idx, val)
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in out), dimension(:) :: arr
    integer, intent(in) :: idx
    integer (INTSIZE), intent(in)  :: val

    integer (INTSIZE), dimension(1) :: val1d

    val1d(1) = val
    call insert (arr, idx, val1d)
end subroutine


subroutine insert_1d_real32 (arr, idx, val, out)
    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:) :: arr
    integer, intent(in) :: idx
    real (PREC), intent(in), dimension(:) :: val
    real (PREC), intent(out), dimension(:) :: out

    include "include/insert_impl.f90"
end subroutine

subroutine insert_inplace_1d_real32 (arr, idx, val)
    integer, parameter :: PREC = real32

    real (PREC), intent(in out), dimension(:) :: arr
    integer, intent(in) :: idx
    real (PREC), intent(in), dimension(:) :: val

    include "include/insert_inplace_impl.f90"
end subroutine

subroutine insert_scalar_real32 (arr, idx, val, out)
    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:) :: arr
    integer, intent(in) :: idx
    real (PREC), intent(in)  :: val
    real (PREC), intent(out), dimension(:) :: out

    real (PREC), dimension(1) :: val1d

    val1d(1) = val
    call insert (arr, idx, val1d, out)
end subroutine

subroutine insert_inplace_scalar_real32 (arr, idx, val)
    integer, parameter :: PREC = real32

    real (PREC), intent(in out), dimension(:) :: arr
    integer, intent(in) :: idx
    real (PREC), intent(in)  :: val

    real (PREC), dimension(1) :: val1d

    val1d(1) = val
    call insert (arr, idx, val1d)
end subroutine


subroutine insert_1d_real64 (arr, idx, val, out)
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:) :: arr
    integer, intent(in) :: idx
    real (PREC), intent(in), dimension(:) :: val
    real (PREC), intent(out), dimension(:) :: out

    include "include/insert_impl.f90"
end subroutine

subroutine insert_inplace_1d_real64 (arr, idx, val)
    integer, parameter :: PREC = real64

    real (PREC), intent(in out), dimension(:) :: arr
    integer, intent(in) :: idx
    real (PREC), intent(in), dimension(:) :: val

    include "include/insert_inplace_impl.f90"
end subroutine

subroutine insert_scalar_real64 (arr, idx, val, out)
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:) :: arr
    integer, intent(in) :: idx
    real (PREC), intent(in)  :: val
    real (PREC), intent(out), dimension(:) :: out

    real (PREC), dimension(1) :: val1d

    val1d(1) = val
    call insert (arr, idx, val1d, out)
end subroutine

subroutine insert_inplace_scalar_real64 (arr, idx, val)
    integer, parameter :: PREC = real64

    real (PREC), intent(in out), dimension(:) :: arr
    integer, intent(in) :: idx
    real (PREC), intent(in)  :: val

    real (PREC), dimension(1) :: val1d

    val1d(1) = val
    call insert (arr, idx, val1d)
end subroutine

end module
