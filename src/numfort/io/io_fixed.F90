

module numfort_io_fixed

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_status
    use numfort_io_common

    implicit none
    private

    public :: read_fixed
    public :: write_fixed

    interface read_fixed
        procedure read_fixed_1d_real32, read_fixed_1d_real64
    end interface

    interface read_fixed
        procedure read_fixed_1d_int32, read_fixed_1d_int64
    end interface

    interface write_fixed
        procedure write_fixed_1d_real32, write_fixed_1d_real64
    end interface

    interface write_fixed
        procedure write_fixed_1d_int32, write_fixed_1d_int64
    end interface

    interface read_fixed
        procedure read_fixed_2d_real32, read_fixed_2d_real64
    end interface

    interface write_fixed
        procedure write_fixed_2d_real32, write_fixed_2d_real64
    end interface

    interface read_fixed
        procedure read_fixed_3d_real32, read_fixed_3d_real64
    end interface

    interface write_fixed
        procedure write_fixed_3d_real32, write_fixed_3d_real64
    end interface

    contains



pure function parse_transform (x) result(res)
    !*  PARSE_TRANSFORM parses a user-provided argument TRANSFORM of character
    !   type and returns the corresponding integer enum.
    character (*), intent(in), optional :: x
    integer (NF_ENUM_KIND) :: res
        !*  Parsed integer value of the transformation that should be applied.
        !   If the argument is not present, the default value is returned.
        !   If the argument is present but contains an invalid value, 0 is
        !   returned.

    character (10) :: lx

    if (.not. present(x)) then
        res = NF_IO_TRANSFORM_NONE
        return
    end if

    lx = x
    call lower (lx)

    select case (lx)
    case ('none')
        res = NF_IO_TRANSFORM_NONE
    case ('transpose')
        res = NF_IO_TRANSFORM_TRANSPOSE
    case ('flatten')
        res = NF_IO_TRANSFORM_FLATTEN
    case default
        res = 0
    end select
end function


!-------------------------------------------------------------------------------
! Routines for 1d-arrays

subroutine read_fixed_1d_real32 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout), dimension(:) :: dat
#include "read_fixed_impl.F90"
end subroutine


subroutine read_fixed_1d_real64 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout), dimension(:) :: dat
#include "read_fixed_impl.F90"
end subroutine


subroutine read_fixed_1d_int32 (path, fmt, dat, transform, status, msg)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(inout), dimension(:) :: dat
#include "read_fixed_impl.F90"
end subroutine


subroutine read_fixed_1d_int64 (path, fmt, dat, transform, status, msg)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(inout), dimension(:) :: dat
#include "read_fixed_impl.F90"
end subroutine


subroutine write_fixed_1d_real32 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:) :: dat
#include "write_fixed_impl.F90"
end subroutine


subroutine write_fixed_1d_real64 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:) :: dat
#include "write_fixed_impl.F90"
end subroutine

subroutine write_fixed_1d_int32 (path, fmt, dat, transform, status, msg)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:) :: dat
#include "write_fixed_impl.F90"
end subroutine


subroutine write_fixed_1d_int64 (path, fmt, dat, transform, status, msg)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:) :: dat
#include "write_fixed_impl.F90"
end subroutine

!-------------------------------------------------------------------------------
! Routines for 2d-arrays

subroutine read_fixed_2d_real32 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout), dimension(:,:) :: dat
#include "read_fixed_2d_impl.F90"
end subroutine


subroutine read_fixed_2d_real64 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout), dimension(:,:) :: dat
#include "read_fixed_2d_impl.F90"
end subroutine


subroutine write_fixed_2d_real32 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout), dimension(:,:) :: dat
#include "write_fixed_2d_impl.F90"
end subroutine


subroutine write_fixed_2d_real64 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout), dimension(:,:) :: dat
#include "write_fixed_2d_impl.F90"
end subroutine


!-------------------------------------------------------------------------------
! Routines for 3d-arrays

subroutine read_fixed_3d_real32 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout), dimension(:,:,:) :: dat
#include "read_fixed_impl.F90"
end subroutine


subroutine read_fixed_3d_real64 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout), dimension(:,:,:) :: dat
#include "read_fixed_impl.F90"
end subroutine


subroutine write_fixed_3d_real32 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(inout), dimension(:,:,:) :: dat
#include "write_fixed_impl.F90"
end subroutine


subroutine write_fixed_3d_real64 (path, fmt, dat, transform, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(inout), dimension(:,:,:) :: dat
#include "write_fixed_impl.F90"
end subroutine

end module
