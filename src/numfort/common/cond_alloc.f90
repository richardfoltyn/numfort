
module numfort_common_cond_alloc
    !*  Module contains implementations of COND_ALLOC routines.

    use, intrinsic :: iso_fortran_env
    use numfort_common_shape, only: has_shape, shape_equal
    implicit none

    private

    public :: cond_alloc

    ! COND_ALLOC for 1d-arrays
    interface cond_alloc
        procedure &
            cond_alloc_1d_int8, cond_alloc_1d_scalar_int8, &
            cond_alloc_1d_scalar_scalar_int8, &
            cond_alloc_1d_int32, cond_alloc_1d_scalar_int32, &
            cond_alloc_1d_scalar_scalar_int32, &
            cond_alloc_1d_int64, cond_alloc_1d_scalar_int64, &
            cond_alloc_1d_scalar_scalar_int64, &
            cond_alloc_1d_real32, cond_alloc_1d_scalar_real32, &
            cond_alloc_1d_scalar_scalar_real32, &
            cond_alloc_1d_real64, cond_alloc_1d_scalar_real64, &
            cond_alloc_1d_scalar_scalar_real64, &
            cond_alloc_1d_logical, cond_alloc_1d_scalar_logical, &
            cond_alloc_1d_scalar_scalar_logical
    end interface

    ! COND_ALLOC for 2d-arrays
    interface cond_alloc
        procedure &
            cond_alloc_2d_int8, cond_alloc_2d_scalar_int8, &
            cond_alloc_2d_int32, cond_alloc_2d_scalar_int32, &
            cond_alloc_2d_int64, cond_alloc_2d_scalar_int64, &
            cond_alloc_2d_real32, cond_alloc_2d_scalar_real32, &
            cond_alloc_2d_real64, cond_alloc_2d_scalar_real64, &
            cond_alloc_2d_logical, cond_alloc_2d_scalar_logical
    end interface


    ! COND_ALLOC for 3d-arrays
    interface cond_alloc
        procedure &
            cond_alloc_3d_int8, cond_alloc_3d_scalar_int8, &
            cond_alloc_3d_int32, cond_alloc_3d_scalar_int32, &
            cond_alloc_3d_int64, cond_alloc_3d_scalar_int64, &
            cond_alloc_3d_real32, cond_alloc_3d_scalar_real32, &
            cond_alloc_3d_real64, cond_alloc_3d_scalar_real64, &
            cond_alloc_3d_logical, cond_alloc_3d_scalar_logical
    end interface


    ! COND_ALLOC for 4d-arrays
    interface cond_alloc
        procedure &
            cond_alloc_4d_int8, cond_alloc_4d_scalar_int8, &
            cond_alloc_4d_int32, cond_alloc_4d_scalar_int32, &
            cond_alloc_4d_int64, cond_alloc_4d_scalar_int64, &
            cond_alloc_4d_real32, cond_alloc_4d_scalar_real32, &
            cond_alloc_4d_real64, cond_alloc_4d_scalar_real64, &
            cond_alloc_4d_logical, cond_alloc_4d_scalar_logical
    end interface

    ! COND_ALLOC for 5d-arrays
    interface cond_alloc
        procedure &
            cond_alloc_5d_int8, cond_alloc_5d_scalar_int8, &
            cond_alloc_5d_int32, cond_alloc_5d_scalar_int32, &
            cond_alloc_5d_int64, cond_alloc_5d_scalar_int64, &
            cond_alloc_5d_real32, cond_alloc_5d_scalar_real32, &
            cond_alloc_5d_real64, cond_alloc_5d_scalar_real64, &
            cond_alloc_5d_logical, cond_alloc_5d_scalar_logical
    end interface
    contains

!-------------------------------------------------------------------------------
! Implementation for 1d-arrays

pure subroutine cond_alloc_1d_scalar_real32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d_scalar.f90"
end subroutine

pure subroutine cond_alloc_1d_scalar_scalar_real32 (arr, shp, source, stat)
    !*  COND_ALLOC wrapper to support scalar shape argument for 1d-arrays.
    integer, parameter :: PREC = real32

    real (PREC), intent(in out), dimension(:), allocatable :: arr
    integer, intent(in) :: shp
    real (PREC), intent(in), optional :: source
    integer, intent(out), optional :: stat

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    call cond_alloc (arr, shp1d, source, stat)
end subroutine

pure subroutine cond_alloc_1d_real32 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d.f90"
end subroutine



pure subroutine cond_alloc_1d_scalar_real64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d_scalar.f90"
end subroutine

pure subroutine cond_alloc_1d_scalar_scalar_real64 (arr, shp, source, stat)
    !*  COND_ALLOC wrapper to support scalar shape argument for 1d-arrays.
    integer, parameter :: PREC = real64

    real (PREC), intent(in out), dimension(:), allocatable :: arr
    integer, intent(in) :: shp
    real (PREC), intent(in), optional :: source
    integer, intent(out), optional :: stat

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    call cond_alloc (arr, shp1d, source, stat)
end subroutine

pure subroutine cond_alloc_1d_real64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d.f90"
end subroutine


pure subroutine cond_alloc_1d_int32 (arr, source, stat)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d.f90"
end subroutine


pure subroutine cond_alloc_1d_scalar_int32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    !
    !   Note: A special implementation is required for INT32 1d-arrays with
    !   all arguments being mandatory, otherwise the interface of several
    !   routines is indistinguishable.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:), allocatable :: arr
    integer, intent(in), dimension(:) :: shp
    integer (INTSIZE), intent(in) :: source
    integer, intent(out) :: stat

    integer :: lstat

    lstat = -1

    if (allocated(arr)) then
        if (.not. has_shape (arr, shp)) deallocate (arr)
    end if

    if (.not. allocated(arr)) then
        allocate (arr(shp(1)), source=source, stat=lstat)
    end if

    stat = lstat
end subroutine

pure subroutine cond_alloc_1d_scalar_scalar_int32 (arr, shp, source, stat)
    !*  COND_ALLOC wrapper to support scalar shape argument for 1d-arrays.
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in out), dimension(:), allocatable :: arr
    integer, intent(in) :: shp
    integer (INTSIZE), intent(in), optional :: source
    integer, intent(out), optional :: stat

    integer :: lstat, lsource
    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    lsource = 0
    if (present(source)) lsource = source
    call cond_alloc (arr, shp1d, lsource, lstat)
    if (present(stat)) stat = lstat
end subroutine


pure subroutine cond_alloc_1d_int64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d.f90"
end subroutine


pure subroutine cond_alloc_1d_scalar_int64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d_scalar.f90"
end subroutine

pure subroutine cond_alloc_1d_scalar_scalar_int64 (arr, shp, source, stat)
    !*  COND_ALLOC wrapper to support scalar shape argument for 1d-arrays.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:), allocatable :: arr
    integer, intent(in) :: shp
    integer (INTSIZE), intent(in), optional :: source
    integer, intent(out), optional :: stat

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    call cond_alloc (arr, shp1d, source, stat)
end subroutine


pure subroutine cond_alloc_1d_int8 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d.f90"
end subroutine


pure subroutine cond_alloc_1d_scalar_int8 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d_scalar.f90"
end subroutine

pure subroutine cond_alloc_1d_scalar_scalar_int8 (arr, shp, source, stat)
    !*  COND_ALLOC wrapper to support scalar shape argument for 1d-arrays.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:), allocatable :: arr
    integer, intent(in) :: shp
    integer (INTSIZE), intent(in), optional :: source
    integer, intent(out), optional :: stat

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    call cond_alloc (arr, shp1d, source, stat)
end subroutine


pure subroutine cond_alloc_1d_scalar_logical (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d_scalar.f90"
end subroutine

pure subroutine cond_alloc_1d_scalar_scalar_logical (arr, shp, source, stat)
    !*  COND_ALLOC wrapper to support scalar shape argument for 1d-arrays.
    logical, intent(in out), dimension(:), allocatable :: arr
    integer, intent(in) :: shp
    logical, intent(in), optional :: source
    integer, intent(out), optional :: stat

    integer, dimension(1) :: shp1d

    shp1d(1) = shp
    call cond_alloc (arr, shp1d, source, stat)
end subroutine

pure subroutine cond_alloc_1d_logical (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), dimension(:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_1d.f90"
end subroutine



!-------------------------------------------------------------------------------
! Implementation for 2d-arrays

pure subroutine cond_alloc_2d_scalar_real32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d_scalar.f90"
end subroutine

pure subroutine cond_alloc_2d_real32 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d.f90"
end subroutine


pure subroutine cond_alloc_2d_scalar_real64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d_scalar.f90"
end subroutine

pure subroutine cond_alloc_2d_real64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d.f90"
end subroutine


pure subroutine cond_alloc_2d_int32 (arr, source, stat)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d.f90"
end subroutine

pure subroutine cond_alloc_2d_scalar_int32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d_scalar.f90"
end subroutine

pure subroutine cond_alloc_2d_int64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d.f90"
end subroutine


pure subroutine cond_alloc_2d_scalar_int64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:,:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d_scalar.f90"
end subroutine


pure subroutine cond_alloc_2d_int8 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d.f90"
end subroutine


pure subroutine cond_alloc_2d_scalar_int8 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:,:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d_scalar.f90"
end subroutine


pure subroutine cond_alloc_2d_scalar_logical (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d_scalar.f90"
end subroutine

pure subroutine cond_alloc_2d_logical (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), dimension(:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_2d.f90"
end subroutine


!-------------------------------------------------------------------------------
! Implementation for 3d-arrays

pure subroutine cond_alloc_3d_scalar_real32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d_scalar.f90"
end subroutine

pure subroutine cond_alloc_3d_real32 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d.f90"
end subroutine


pure subroutine cond_alloc_3d_scalar_real64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d_scalar.f90"
end subroutine

pure subroutine cond_alloc_3d_real64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d.f90"
end subroutine


pure subroutine cond_alloc_3d_int32 (arr, source, stat)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d.f90"
end subroutine

pure subroutine cond_alloc_3d_scalar_int32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d_scalar.f90"
end subroutine

pure subroutine cond_alloc_3d_int64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d.f90"
end subroutine


pure subroutine cond_alloc_3d_scalar_int64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:,:,:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d_scalar.f90"
end subroutine


pure subroutine cond_alloc_3d_int8 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d.f90"
end subroutine


pure subroutine cond_alloc_3d_scalar_int8 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:,:,:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d_scalar.f90"
end subroutine


pure subroutine cond_alloc_3d_scalar_logical (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d_scalar.f90"
end subroutine

pure subroutine cond_alloc_3d_logical (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), dimension(:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_3d.f90"
end subroutine


!-------------------------------------------------------------------------------
! Implementation for 4d-arrays

pure subroutine cond_alloc_4d_scalar_real32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d_scalar.f90"
end subroutine

pure subroutine cond_alloc_4d_real32 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d.f90"
end subroutine


pure subroutine cond_alloc_4d_scalar_real64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d_scalar.f90"
end subroutine

pure subroutine cond_alloc_4d_real64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d.f90"
end subroutine


pure subroutine cond_alloc_4d_int32 (arr, source, stat)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d.f90"
end subroutine

pure subroutine cond_alloc_4d_scalar_int32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d_scalar.f90"
end subroutine

pure subroutine cond_alloc_4d_int64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d.f90"
end subroutine


pure subroutine cond_alloc_4d_scalar_int64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:,:,:,:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d_scalar.f90"
end subroutine


pure subroutine cond_alloc_4d_int8 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d.f90"
end subroutine


pure subroutine cond_alloc_4d_scalar_int8 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:,:,:,:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d_scalar.f90"
end subroutine


pure subroutine cond_alloc_4d_scalar_logical (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d_scalar.f90"
end subroutine

pure subroutine cond_alloc_4d_logical (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), dimension(:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_4d.f90"
end subroutine

!-------------------------------------------------------------------------------
! Implementation for 5d-arrays

pure subroutine cond_alloc_5d_scalar_real32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d_scalar.f90"
end subroutine

pure subroutine cond_alloc_5d_real32 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real32
    real (PREC), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:,:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d.f90"
end subroutine


pure subroutine cond_alloc_5d_scalar_real64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d_scalar.f90"
end subroutine

pure subroutine cond_alloc_5d_real64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: PREC = real64
    real (PREC), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    real (PREC), intent(in), dimension(:,:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d.f90"
end subroutine


pure subroutine cond_alloc_5d_int32 (arr, source, stat)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d.f90"
end subroutine

pure subroutine cond_alloc_5d_scalar_int32 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d_scalar.f90"
end subroutine

pure subroutine cond_alloc_5d_int64 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d.f90"
end subroutine


pure subroutine cond_alloc_5d_scalar_int64 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d_scalar.f90"
end subroutine


pure subroutine cond_alloc_5d_int8 (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d.f90"
end subroutine


pure subroutine cond_alloc_5d_scalar_int8 (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in out), dimension(:,:,:,:,:), allocatable :: arr
    integer (INTSIZE), intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d_scalar.f90"
end subroutine


pure subroutine cond_alloc_5d_scalar_logical (arr, shp, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), optional :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d_scalar.f90"
end subroutine

pure subroutine cond_alloc_5d_logical (arr, source, stat)
    !*  COND_ALLOC conditinally allocates a given array if it is not already
    !   allocated, or its allocated shape differs from the desired shape
    !   or source array.
    logical, intent(in out), dimension(:,:,:,:,:), allocatable :: arr
        !*  Array to be conditionally allocated
    logical, intent(in), dimension(:,:,:,:,:) :: source
        !*  If present, used as source value if array ARR needs to be
        !   (re)allocated.
    include "include/cond_alloc_5d.f90"
end subroutine


end module
