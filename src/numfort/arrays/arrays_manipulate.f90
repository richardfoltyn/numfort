

module numfort_arrays_manipulate
    !*  Module contains various array manipulation routines that do not
    !   fit any more specific category.

    use, intrinsic :: iso_fortran_env
    use numfort_common_status
    use numfort_arrays_create, only: arange
    use numfort_arrays_setops, only: unique, setdiff

    implicit none
    private

    public :: insert
    public :: delete


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

    ! DELETE for 1d-arrays
    interface delete
        procedure &
            delete_1d_real32, delete_1d_scalar_real32, &
            delete_1d_real64, delete_1d_scalar_real64, &
            delete_1d_int32, delete_1d_scalar_int32, &
            delete_1d_int64, delete_1d_scalar_int64
    end interface

    ! DELETE for 2d-arrays
    interface delete
        procedure &
            delete_2d_real32, delete_2d_scalar_real32, &
            delete_2d_real64, delete_2d_scalar_real64, &
            delete_2d_int32, delete_2d_scalar_int32, &
            delete_2d_int64, delete_2d_scalar_int64
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

!-------------------------------------------------------------------------------
! DELETE routines for 1d-arrays

subroutine delete_1d_real32 (arr, idx, out, dim, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the elements at
    !   positions specified by IDX.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:) :: arr
        !*  Input array.
    real (PREC), intent(out), dimension(:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT is
        !   not large enough to hold all non-deleted elements from ARR, a
        !   truncated subset of those elements is returned.
    include "include/delete_1d_impl.f90"
end subroutine

subroutine delete_1d_scalar_real32 (arr, idx, out, dim, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the element at
    !   position specified by IDX.
    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:) :: arr
        !*  Input array.
    integer, intent(in) :: idx
        !*  Element which should be deleted from ARR.
    real (PREC), intent(out), dimension(:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT is
        !   not large enough to hold all non-deleted elements from ARR, a
        !   truncated subset of those elements is returned.
    integer, intent(in), optional :: dim
        !*  The dimension along with elements should be deleted. Ignored for
        !   1d-arrays. An error is raised if this argument is present and DIM /= 1.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT that holds
        !   valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer, dimension(1) :: idx1d
    idx1d(1) = idx
    call delete (arr, idx1d, out, dim, .true., n, status)
end subroutine


subroutine delete_1d_real64 (arr, idx, out, dim, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the elements at
    !   positions specified by IDX.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:) :: arr
        !*  Input array.
    real (PREC), intent(out), dimension(:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT is
        !   not large enough to hold all non-deleted elements from ARR, a
        !   truncated subset of those elements is returned.
    include "include/delete_1d_impl.f90"
end subroutine

subroutine delete_1d_scalar_real64 (arr, idx, out, dim, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the element at
    !   position specified by IDX.
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:) :: arr
        !*  Input array.
    integer, intent(in) :: idx
        !*  Element which should be deleted from ARR.
    real (PREC), intent(out), dimension(:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT is
        !   not large enough to hold all non-deleted elements from ARR, a
        !   truncated subset of those elements is returned.
    integer, intent(in), optional :: dim
        !*  The dimension along with elements should be deleted. Ignored for
        !   1d-arrays. An error is raised if this argument is present and DIM /= 1.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT that holds
        !   valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer, dimension(1) :: idx1d
    idx1d(1) = idx
    call delete (arr, idx1d, out, dim, .true., n, status)
end subroutine


subroutine delete_1d_int32 (arr, idx, out, dim, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the elements at
    !   positions specified by IDX.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array.
    integer (INTSIZE), intent(out), dimension(:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT is
        !   not large enough to hold all non-deleted elements from ARR, a
        !   truncated subset of those elements is returned.
    include "include/delete_1d_impl.f90"
end subroutine

subroutine delete_1d_scalar_int32 (arr, idx, out, dim, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the element at
    !   position specified by IDX.
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array.
    integer, intent(in) :: idx
        !*  Element which should be deleted from ARR.
    integer (INTSIZE), intent(out), dimension(:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT is
        !   not large enough to hold all non-deleted elements from ARR, a
        !   truncated subset of those elements is returned.
    integer, intent(in), optional :: dim
        !*  The dimension along with elements should be deleted. Ignored for
        !   1d-arrays. An error is raised if this argument is present and DIM /= 1.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT that holds
        !   valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer, dimension(1) :: idx1d
    idx1d(1) = idx
    call delete (arr, idx1d, out, dim, .true., n, status)
end subroutine


subroutine delete_1d_int64 (arr, idx, out, dim, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the elements at
    !   positions specified by IDX.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array.
    integer (INTSIZE), intent(out), dimension(:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT is
        !   not large enough to hold all non-deleted elements from ARR, a
        !   truncated subset of those elements is returned.
    include "include/delete_1d_impl.f90"
end subroutine

subroutine delete_1d_scalar_int64 (arr, idx, out, dim, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the element at
    !   position specified by IDX.
    integer, parameter :: INTSIZE = int64

    integer (INTSIZE), intent(in), dimension(:) :: arr
        !*  Input array.
    integer, intent(in) :: idx
        !*  Element which should be deleted from ARR.
    integer (INTSIZE), intent(out), dimension(:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT is
        !   not large enough to hold all non-deleted elements from ARR, a
        !   truncated subset of those elements is returned.
    integer, intent(in), optional :: dim
        !*  The dimension along with elements should be deleted. Ignored for
        !   1d-arrays. An error is raised if this argument is present and DIM /= 1.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT that holds
        !   valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer, dimension(1) :: idx1d
    idx1d(1) = idx
    call delete (arr, idx1d, out, dim, .true., n, status)
end subroutine

!-------------------------------------------------------------------------------
! DELETE routines for 2d-arrays

subroutine delete_2d_real32 (arr, idx, dim, out, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the elements at
    !   positions specified by IDX along the dimension DIM.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:) :: arr
        !*  Input array.
    real (PREC), intent(out), dimension(:,:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT
        !   is not large enough to hold all non-deleted elements of ARR,
        !   the routine exits, leaving OUT unmodified.
    include "include/delete_2d_impl.f90"
end subroutine

subroutine delete_2d_scalar_real32 (arr, idx, dim, out, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the element at
    !   position specified by IDX along the dimension DIM.
    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:,:) :: arr
        !*  Input array.
    integer, intent(in) :: idx
        !*  Element which should be deleted from ARR.
    integer, intent(in) :: dim
        !*  The dimension along with elements should be deleted.
    real (PREC), intent(out), dimension(:,:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT
        !   is not large enough to hold all non-deleted elements of ARR,
        !   the routine exits, leaving OUT unmodified.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT, dimension DIM,
        !   that holds valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer, dimension(1) :: idx1d
    idx1d(1) = idx
    call delete (arr, idx1d, dim, out, .true., n, status)
end subroutine


subroutine delete_2d_real64 (arr, idx, dim, out, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the elements at
    !   positions specified by IDX along the dimension DIM.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:) :: arr
        !*  Input array.
    real (PREC), intent(out), dimension(:,:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT
        !   is not large enough to hold all non-deleted elements of ARR,
        !   the routine exits, leaving OUT unmodified.
    include "include/delete_2d_impl.f90"
end subroutine

subroutine delete_2d_scalar_real64 (arr, idx, dim, out, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the element at
    !   position specified by IDX along the dimension DIM.
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:,:) :: arr
        !*  Input array.
    integer, intent(in) :: idx
        !*  Element which should be deleted from ARR.
    integer, intent(in) :: dim
        !*  The dimension along with elements should be deleted.
    real (PREC), intent(out), dimension(:,:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT
        !   is not large enough to hold all non-deleted elements of ARR,
        !   the routine exits, leaving OUT unmodified.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT, dimension DIM,
        !   that holds valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer, dimension(1) :: idx1d
    idx1d(1) = idx
    call delete (arr, idx1d, dim, out, .true., n, status)
end subroutine


subroutine delete_2d_int32 (arr, idx, dim, out, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the elements at
    !   positions specified by IDX along the dimension DIM.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:) :: arr
        !*  Input array.
    integer (INTSIZE), intent(out), dimension(:,:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT
        !   is not large enough to hold all non-deleted elements of ARR,
        !   the routine exits, leaving OUT unmodified.
    include "include/delete_2d_impl.f90"
end subroutine

subroutine delete_2d_scalar_int32 (arr, idx, dim, out, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the element at
    !   position specified by IDX along the dimension DIM.
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:,:) :: arr
        !*  Input array.
    integer, intent(in) :: idx
        !*  Element which should be deleted from ARR.
    integer, intent(in) :: dim
        !*  The dimension along with elements should be deleted.
    integer (INTSIZE), intent(out), dimension(:,:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT
        !   is not large enough to hold all non-deleted elements of ARR,
        !   the routine exits, leaving OUT unmodified.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT, dimension DIM,
        !   that holds valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer, dimension(1) :: idx1d
    idx1d(1) = idx
    call delete (arr, idx1d, dim, out, .true., n, status)
end subroutine


subroutine delete_2d_int64 (arr, idx, dim, out, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the elements at
    !   positions specified by IDX along the dimension DIM.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:) :: arr
        !*  Input array.
    integer (INTSIZE), intent(out), dimension(:,:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT
        !   is not large enough to hold all non-deleted elements of ARR,
        !   the routine exits, leaving OUT unmodified.
    include "include/delete_2d_impl.f90"
end subroutine

subroutine delete_2d_scalar_int64 (arr, idx, dim, out, sorted, n, status)
    !*  DELETE returns the contents of array ARR without the element at
    !   position specified by IDX along the dimension DIM.
    integer, parameter :: INTSIZE = int64

    integer (INTSIZE), intent(in), dimension(:,:) :: arr
        !*  Input array.
    integer, intent(in) :: idx
        !*  Element which should be deleted from ARR.
    integer, intent(in) :: dim
        !*  The dimension along with elements should be deleted.
    integer (INTSIZE), intent(out), dimension(:,:) :: out
        !*  A copy of ARR with elements specified by IDX removed. If OUT
        !   is not large enough to hold all non-deleted elements of ARR,
        !   the routine exits, leaving OUT unmodified.
    logical, intent(in), optional :: sorted
        !*  If present and true, the values in IDX are assumed to be unique and
        !   sorted in ascending order. Otherwise elements will be sorted
        !   by the routine.
    integer, intent(out), optional :: n
        !*  If present, contains the highest index on array OUT, dimension DIM,
        !   that holds valid data.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status code.

    integer, dimension(1) :: idx1d
    idx1d(1) = idx
    call delete (arr, idx1d, dim, out, .true., n, status)
end subroutine


end module
