
module numfort_arrays_shape

    use, intrinsic :: iso_fortran_env
    
    implicit none
    
    private 
    
    public :: shape_equal
    
    interface shape_equal
        module procedure shape_equal_1d_real32, shape_equal_1d_real64, &
            shape_equal_1d_int32, shape_equal_1d_int64
    end interface
    
    ! 2d-array routines
    interface shape_equal
        module procedure shape_equal_2d_real32, shape_equal_2d_real64, &
            shape_equal_2d_int32, shape_equal_2d_int64
    end interface
    
    ! 3d-array routines
    interface shape_equal
        module procedure shape_equal_3d_real32, shape_equal_3d_real64, &
            shape_equal_3d_int32, shape_equal_3d_int64
    end interface
    
    ! 4d-array routines
    interface shape_equal
        module procedure shape_equal_4d_real32, shape_equal_4d_real64, &
            shape_equal_4d_int32, shape_equal_4d_int64
    end interface
contains

pure function shape_equal_1d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: ND = 1
    real (PREC), intent(in), dimension(:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_1d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: ND = 1
    real (PREC), intent(in), dimension(:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_1d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: ND = 1
    integer (INTSIZE), intent(in), dimension(:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_1d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: ND = 1
    integer (INTSIZE), intent(in), dimension(:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_2d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: ND = 2
    real (PREC), intent(in), dimension(:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_2d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: ND = 2
    real (PREC), intent(in), dimension(:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_2d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: ND = 2
    integer (INTSIZE), intent(in), dimension(:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_2d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: ND = 2
    integer (INTSIZE), intent(in), dimension(:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_3d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: ND = 3
    real (PREC), intent(in), dimension(:,:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_3d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: ND = 3
    real (PREC), intent(in), dimension(:,:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_3d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: ND = 3
    integer (INTSIZE), intent(in), dimension(:,:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_3d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: ND = 3
    integer (INTSIZE), intent(in), dimension(:,:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_4d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: ND = 4
    real (PREC), intent(in), dimension(:,:,:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_4d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: ND = 4
    real (PREC), intent(in), dimension(:,:,:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_4d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: ND = 4
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_4d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: ND = 4
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: arr1, arr2
    
    include "include/shape_equal_impl.f90"
end function

end module
