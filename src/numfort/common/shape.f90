
module numfort_common_shape

    use, intrinsic :: iso_fortran_env

    implicit none

    private

    public :: shape_equal
    public :: has_shape

    interface shape_equal
        procedure &
            shape_equal_1d_real32, &
            shape_equal_1d_real64, &
            shape_equal_1d_int8, &
            shape_equal_1d_int32, &
            shape_equal_1d_int64, &
            shape_equal_1d_logical
    end interface

    ! 2d-array routines
    interface shape_equal
        procedure &
            shape_equal_2d_real32, &
            shape_equal_2d_real64, &
            shape_equal_2d_int8, &
            shape_equal_2d_int32, &
            shape_equal_2d_int64, &
            shape_equal_2d_logical
    end interface

    ! 3d-array routines
    interface shape_equal
        procedure &
            shape_equal_3d_real32, &
            shape_equal_3d_real64, &
            shape_equal_3d_int8, &
            shape_equal_3d_int32, &
            shape_equal_3d_int64, &
            shape_equal_3d_logical
    end interface

    ! 4d-array routines
    interface shape_equal
        procedure &
            shape_equal_4d_real32, &
            shape_equal_4d_real64, &
            shape_equal_4d_int8, &
            shape_equal_4d_int32, &
            shape_equal_4d_int64, &
            shape_equal_4d_logical
    end interface

    ! 5d-array routines
    interface shape_equal
        procedure &
            shape_equal_5d_real32, &
            shape_equal_5d_real64, &
            shape_equal_5d_int8, &
            shape_equal_5d_int32, &
            shape_equal_5d_int64, &
            shape_equal_5d_logical
    end interface

    ! 6d-array routines
    interface shape_equal
        procedure &
            shape_equal_6d_real32, &
            shape_equal_6d_real64, &
            shape_equal_6d_int8, &
            shape_equal_6d_int32, &
            shape_equal_6d_int64, &
            shape_equal_6d_logical
    end interface

    ! 1d-array routines
    interface has_shape
        procedure &
            has_shape_1d_logical, has_shape_1d_scalar_logical, &
            has_shape_1d_int8, has_shape_1d_scalar_int8, &
            has_shape_1d_int32, has_shape_1d_scalar_int32, &
            has_shape_1d_int64, has_shape_1d_scalar_int64, &
            has_shape_1d_real32, has_shape_1d_scalar_real32, &
            has_shape_1d_real64, has_shape_1d_scalar_real64
    end interface

    ! 2d-array routines
    interface has_shape
        procedure &
            has_shape_2d_logical, &
            has_shape_2d_int8, &
            has_shape_2d_int32, &
            has_shape_2d_int64, &
            has_shape_2d_real32, &
            has_shape_2d_real64
    end interface

    ! 3d-array routines
    interface has_shape
        procedure &
            has_shape_3d_logical, &
            has_shape_3d_int8, &
            has_shape_3d_int32, &
            has_shape_3d_int64, &
            has_shape_3d_real32, &
            has_shape_3d_real64
    end interface

    ! 4d-array routines
    interface has_shape
        procedure &
            has_shape_4d_logical, &
            has_shape_4d_int8, &
            has_shape_4d_int32, &
            has_shape_4d_int64, &
            has_shape_4d_real32, &
            has_shape_4d_real64
    end interface

    ! 5d-array routines
    interface has_shape
        procedure &
            has_shape_5d_logical, &
            has_shape_5d_int8, &
            has_shape_5d_int32, &
            has_shape_5d_int64, &
            has_shape_5d_real32, &
            has_shape_5d_real64
    end interface

    ! 6d-array routines
    interface has_shape
        procedure &
            has_shape_6d_logical, &
            has_shape_6d_int8, &
            has_shape_6d_int32, &
            has_shape_6d_int64, &
            has_shape_6d_real32, &
            has_shape_6d_real64
    end interface

    contains

!-------------------------------------------------------------------------------
! SHAPE_EQUAL for 1d-arrays

pure function shape_equal_1d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: NDIM = 1
    real (PREC), intent(in), dimension(:) :: arr1
    real (PREC), intent(in), dimension(:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_1d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: NDIM = 1
    real (PREC), intent(in), dimension(:) :: arr1
    real (PREC), intent(in), dimension(:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_1d_int8 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 1
    integer (INTSIZE), intent(in), dimension(:) :: arr1
    integer (INTSIZE), intent(in), dimension(:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_1d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: NDIM = 1
    integer (INTSIZE), intent(in), dimension(:) :: arr1
    integer (INTSIZE), intent(in), dimension(:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_1d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: NDIM = 1
    integer (INTSIZE), intent(in), dimension(:) :: arr1
    integer (INTSIZE), intent(in), dimension(:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_1d_logical (arr1, arr2) result(res)
    integer, parameter :: NDIM = 1
    logical, intent(in), dimension(:) :: arr1
    logical, intent(in), dimension(:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

!-------------------------------------------------------------------------------
! SHAPE_EQUAL for 2d-arrays

pure function shape_equal_2d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: NDIM = 2
    real (PREC), intent(in), dimension(:,:) :: arr1
    real (PREC), intent(in), dimension(:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_2d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: NDIM = 2
    real (PREC), intent(in), dimension(:,:) :: arr1
    real (PREC), intent(in), dimension(:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_2d_int8 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 2
    integer (INTSIZE), intent(in), dimension(:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_2d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: NDIM = 2
    integer (INTSIZE), intent(in), dimension(:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_2d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: NDIM = 2
    integer (INTSIZE), intent(in), dimension(:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_2d_logical (arr1, arr2) result(res)
    integer, parameter :: NDIM = 2
    logical, intent(in), dimension(:,:) :: arr1
    logical, intent(in), dimension(:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

!-------------------------------------------------------------------------------
! SHAPE_EQUAL for 3d-arrays

pure function shape_equal_3d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: NDIM = 3
    real (PREC), intent(in), dimension(:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_3d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: NDIM = 3
    real (PREC), intent(in), dimension(:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_3d_int8 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 3
    integer (INTSIZE), intent(in), dimension(:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_3d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: NDIM = 3
    integer (INTSIZE), intent(in), dimension(:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_3d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: NDIM = 3
    integer (INTSIZE), intent(in), dimension(:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_3d_logical (arr1, arr2) result(res)
    integer, parameter :: NDIM = 3
    logical, intent(in), dimension(:,:,:) :: arr1
    logical, intent(in), dimension(:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

!-------------------------------------------------------------------------------
! SHAPE_EQUAL for 4d-arrays

pure function shape_equal_4d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: NDIM = 4
    real (PREC), intent(in), dimension(:,:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_4d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: NDIM = 4
    real (PREC), intent(in), dimension(:,:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_4d_int8 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 4
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_4d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: NDIM = 4
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_4d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: NDIM = 4
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_4d_logical (arr1, arr2) result(res)
    integer, parameter :: NDIM = 4
    logical, intent(in), dimension(:,:,:,:) :: arr1
    logical, intent(in), dimension(:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function


!-------------------------------------------------------------------------------
! SHAPE_EQUAL for 5d-arrays

pure function shape_equal_5d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: NDIM = 5
    real (PREC), intent(in), dimension(:,:,:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_5d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: NDIM = 5
    real (PREC), intent(in), dimension(:,:,:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_5d_int8 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 5
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_5d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: NDIM = 5
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_5d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: NDIM = 5
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_5d_logical (arr1, arr2) result(res)
    integer, parameter :: NDIM = 5
    logical, intent(in), dimension(:,:,:,:,:) :: arr1
    logical, intent(in), dimension(:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function


!-------------------------------------------------------------------------------
! SHAPE_EQUAL for 6d-arrays

pure function shape_equal_6d_real32 (arr1, arr2) result(res)
    integer, parameter :: PREC = real32
    integer, parameter :: NDIM = 6
    real (PREC), intent(in), dimension(:,:,:,:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_6d_real64 (arr1, arr2) result(res)
    integer, parameter :: PREC = real64
    integer, parameter :: NDIM = 6
    real (PREC), intent(in), dimension(:,:,:,:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_6d_int8 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 6
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_6d_int32 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: NDIM = 6
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_6d_int64 (arr1, arr2) result(res)
    integer, parameter :: INTSIZE = int64
    integer, parameter :: NDIM = 6
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:) :: arr1
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function

pure function shape_equal_6d_logical (arr1, arr2) result(res)
    integer, parameter :: NDIM = 6
    logical, intent(in), dimension(:,:,:,:,:,:) :: arr1
    logical, intent(in), dimension(:,:,:,:,:,:), optional :: arr2

    include "include/shape_equal_impl.f90"
end function



!-------------------------------------------------------------------------------
! HAS_SHAPE for 1d-arrays

pure function has_shape_1d_real32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 1
    include "include/has_shape_impl.f90"
end function

pure function has_shape_1d_scalar_real32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape
    !   (convenience wrapper for 1d-arrays that permits scalar SHP arguments).
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, intent(in) :: shp
        !*  Shape which ARR shape should be compared to. Passed as scalar
        !   for 1d arrays.
    logical :: res
        !*  True if given array has the desired shape and false otherwise.

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    res = has_shape (arr, shp1d)
end function

pure function has_shape_1d_real64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 1
    include "include/has_shape_impl.f90"
end function

pure function has_shape_1d_scalar_real64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape
    !   (convenience wrapper for 1d-arrays that permits scalar SHP arguments).
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, intent(in) :: shp
        !*  Shape which ARR shape should be compared to. Passed as scalar
        !   for 1d arrays.
    logical :: res
        !*  True if given array has the desired shape and false otherwise.

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    res = has_shape (arr, shp1d)
end function

pure function has_shape_1d_int8 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 1
    include "include/has_shape_impl.f90"
end function

pure function has_shape_1d_scalar_int8 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape
    !   (convenience wrapper for 1d-arrays that permits scalar SHP arguments).
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, intent(in) :: shp
        !*  Shape which ARR shape should be compared to. Passed as scalar
        !   for 1d arrays.
    logical :: res
        !*  True if given array has the desired shape and false otherwise.

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    res = has_shape (arr, shp1d)
end function

pure function has_shape_1d_int32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 1
    include "include/has_shape_impl.f90"
end function

pure function has_shape_1d_scalar_int32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape
    !   (convenience wrapper for 1d-arrays that permits scalar SHP arguments).
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, intent(in) :: shp
        !*  Shape which ARR shape should be compared to. Passed as scalar
        !   for 1d arrays.
    logical :: res
        !*  True if given array has the desired shape and false otherwise.

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    res = has_shape (arr, shp1d)
end function

pure function has_shape_1d_int64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 1
    include "include/has_shape_impl.f90"
end function

pure function has_shape_1d_scalar_int64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape
    !   (convenience wrapper for 1d-arrays that permits scalar SHP arguments).
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, intent(in) :: shp
        !*  Shape which ARR shape should be compared to. Passed as scalar
        !   for 1d arrays.
    logical :: res
        !*  True if given array has the desired shape and false otherwise.

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    res = has_shape (arr, shp1d)
end function

pure function has_shape_1d_logical (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    logical, intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 1
    include "include/has_shape_impl.f90"
end function

pure function has_shape_1d_scalar_logical (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape
    !   (convenience wrapper for 1d-arrays that permits scalar SHP arguments).
    logical, intent(in), dimension(:) :: arr
        !*  Input array to check
    integer, intent(in) :: shp
        !*  Shape which ARR shape should be compared to. Passed as scalar
        !   for 1d arrays.
    logical :: res
        !*  True if given array has the desired shape and false otherwise.

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    res = has_shape (arr, shp1d)
end function

!-------------------------------------------------------------------------------
! HAS_SHAPE for 2d-arrays

pure function has_shape_2d_real32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 2
    include "include/has_shape_impl.f90"
end function

pure function has_shape_2d_real64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 2
    include "include/has_shape_impl.f90"
end function

pure function has_shape_2d_int8 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 2
    include "include/has_shape_impl.f90"
end function

pure function has_shape_2d_int32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 2
    include "include/has_shape_impl.f90"
end function

pure function has_shape_2d_int64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 2
    include "include/has_shape_impl.f90"
end function

pure function has_shape_2d_logical (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    logical, intent(in), dimension(:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 2
    include "include/has_shape_impl.f90"
end function

!-------------------------------------------------------------------------------
! HAS_SHAPE for 3d-arrays

pure function has_shape_3d_real32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 3
    include "include/has_shape_impl.f90"
end function

pure function has_shape_3d_real64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 3
    include "include/has_shape_impl.f90"
end function

pure function has_shape_3d_int8 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 3
    include "include/has_shape_impl.f90"
end function

pure function has_shape_3d_int32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 3
    include "include/has_shape_impl.f90"
end function

pure function has_shape_3d_int64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 3
    include "include/has_shape_impl.f90"
end function

pure function has_shape_3d_logical (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    logical, intent(in), dimension(:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 3
    include "include/has_shape_impl.f90"
end function


!-------------------------------------------------------------------------------
! HAS_SHAPE for 4d-arrays

pure function has_shape_4d_real32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 4
    include "include/has_shape_impl.f90"
end function

pure function has_shape_4d_real64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 4
    include "include/has_shape_impl.f90"
end function

pure function has_shape_4d_int8 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 4
    include "include/has_shape_impl.f90"
end function

pure function has_shape_4d_int32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 4
    include "include/has_shape_impl.f90"
end function

pure function has_shape_4d_int64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 4
    include "include/has_shape_impl.f90"
end function

pure function has_shape_4d_logical (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    logical, intent(in), dimension(:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 4
    include "include/has_shape_impl.f90"
end function

!-------------------------------------------------------------------------------
! HAS_SHAPE for 5d-arrays

pure function has_shape_5d_real32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 5
    include "include/has_shape_impl.f90"
end function

pure function has_shape_5d_real64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 5
    include "include/has_shape_impl.f90"
end function

pure function has_shape_5d_int8 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 5
    include "include/has_shape_impl.f90"
end function

pure function has_shape_5d_int32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 5
    include "include/has_shape_impl.f90"
end function

pure function has_shape_5d_int64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 5
    include "include/has_shape_impl.f90"
end function

pure function has_shape_5d_logical (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    logical, intent(in), dimension(:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 5
    include "include/has_shape_impl.f90"
end function

!-------------------------------------------------------------------------------
! HAS_SHAPE for 6d-arrays

pure function has_shape_6d_real32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 6
    include "include/has_shape_impl.f90"
end function

pure function has_shape_6d_real64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 6
    include "include/has_shape_impl.f90"
end function

pure function has_shape_6d_int8 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 6
    include "include/has_shape_impl.f90"
end function

pure function has_shape_6d_int32 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 6
    include "include/has_shape_impl.f90"
end function

pure function has_shape_6d_int64 (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 6
    include "include/has_shape_impl.f90"
end function

pure function has_shape_6d_logical (arr, shp) result(res)
    !*  HAS_SHAPE verifies that an array argument has a given shape.
    logical, intent(in), dimension(:,:,:,:,:,:) :: arr
        !*  Input array to check
    integer, parameter :: NDIM = 6
    include "include/has_shape_impl.f90"
end function


end module
