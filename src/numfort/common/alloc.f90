
module numfort_common_alloc
    !*  Module with helper routines to help management of dynamically allocated
    !   objects.

    use, intrinsic :: iso_fortran_env
    use numfort_common_shape, only: has_shape, shape_equal
    implicit none

    private

    public :: assert_alloc_ptr, assert_dealloc_ptr
    public :: cond_alloc

    interface assert_alloc_ptr
        module procedure assert_alloc_ptr_1d_real32, assert_alloc_ptr_1d_real64
    end interface

    interface assert_dealloc_ptr
        module procedure assert_dealloc_ptr_1d_real32, assert_dealloc_ptr_1d_real64
    end interface


    ! COND_ALLOC for 1d-arrays
    interface cond_alloc
        procedure cond_alloc_1d_int32, cond_alloc_1d_scalar_int32, &
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
        procedure cond_alloc_2d_int32, cond_alloc_2d_scalar_int32, &
            cond_alloc_2d_int64, cond_alloc_2d_scalar_int64, &
            cond_alloc_2d_real32, cond_alloc_2d_scalar_real32, &
            cond_alloc_2d_real64, cond_alloc_2d_scalar_real64, &
            cond_alloc_2d_logical, cond_alloc_2d_scalar_logical
    end interface


    ! COND_ALLOC for 3d-arrays
    interface cond_alloc
        procedure cond_alloc_3d_int32, cond_alloc_3d_scalar_int32, &
            cond_alloc_3d_int64, cond_alloc_3d_scalar_int64, &
            cond_alloc_3d_real32, cond_alloc_3d_scalar_real32, &
            cond_alloc_3d_real64, cond_alloc_3d_scalar_real64, &
            cond_alloc_3d_logical, cond_alloc_3d_scalar_logical
    end interface


    ! COND_ALLOC for 4d-arrays
    interface cond_alloc
        procedure cond_alloc_4d_int32, cond_alloc_4d_scalar_int32, &
            cond_alloc_4d_int64, cond_alloc_4d_scalar_int64, &
            cond_alloc_4d_real32, cond_alloc_4d_scalar_real32, &
            cond_alloc_4d_real64, cond_alloc_4d_scalar_real64, &
            cond_alloc_4d_logical, cond_alloc_4d_scalar_logical
    end interface

    contains

! ------------------------------------------------------------------------------
! ASSERT_ALLOC_PTR routines

subroutine assert_alloc_ptr_1d_real32 (x, n, ptr_x)
    !*  ASSERT_ALLOC_PTR ensures that ptr_x points to allocated memory:
    !   If argument x is present, ptr_x points to x and no additional
    !   memory is allocated. If x is not present, then ptr_x points to a newly
    !   allocated block of memory of size n.

    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:), target, optional :: x
    integer, intent(in) :: n
    real (PREC), intent(out), dimension(:), pointer :: ptr_x

    if (present(x)) then
        if (size(x) >= n) then
            ptr_x => x
        else
            allocate (ptr_x(n))
        end if
    else
        allocate (ptr_x(n))
    end if
end subroutine

subroutine assert_alloc_ptr_1d_real64 (x, n, ptr_x)
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:), target, optional :: x
    integer, intent(in) :: n
    real (PREC), intent(out), dimension(:), pointer :: ptr_x

    if (present(x)) then
        if (size(x) >= n) then
            ptr_x => x
        else
            allocate (ptr_x(n))
        end if
    else
        allocate (ptr_x(n))
    end if
end subroutine

! ------------------------------------------------------------------------------
! ASSERT_DEALLOC_PTR routines

pure subroutine assert_dealloc_ptr_1d_real32 (x, ptr_x)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_x
    !   needs to be released only if argument x is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_x.

    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:), target, optional :: x
    real (PREC), intent(in out), dimension(:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine

pure subroutine assert_dealloc_ptr_1d_real64 (x, ptr_x)
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:), target, optional :: x
    real (PREC), intent(in out), dimension(:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine

!-------------------------------------------------------------------------------
! COND_ALLOC routines

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


end module
