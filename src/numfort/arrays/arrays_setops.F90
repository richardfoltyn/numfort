module numfort_arrays_setops

    use, intrinsic :: iso_fortran_env
    use m_unirnk, only: unirnk

    use numfort_common_status
    use numfort_common_input_checks
    use numfort_arrays_sort

    implicit none

    private

    public :: unique, setdiff, intersect, union
    public :: split_uniform
    public :: operator(.in.)

    interface unique
        procedure unique_real32, unique_real64, unique_int32
    end interface

    interface setdiff
        procedure setdiff_real32, setdiff_real64, setdiff_int32
    end interface

    interface intersect
        procedure intersect_real32, intersect_real64, intersect_int32
    end interface

    interface union
        procedure union_real32, union_real64, union_int32
    end interface

    interface split_uniform
        procedure split_uniform_int32
    end interface

    interface operator(.in.)
        procedure set_in_int32, set_in_int64
    end interface

    contains

! ------------------------------------------------------------------------------
! UNIQUE routines

subroutine unique_real32 (arr, n, uarr, idx)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: arr
    real (PREC), intent(out), dimension(:), optional :: uarr

#include "unique_impl.f90"
end subroutine


subroutine unique_real64 (arr, n, uarr, idx)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: arr
    real (PREC), intent(out), dimension(:), optional :: uarr

#include "unique_impl.f90"
end subroutine


subroutine unique_int32 (arr, n, uarr, idx)
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:) :: arr
    integer (INTSIZE), intent(out), dimension(:), optional :: uarr

#include "unique_impl.f90"
end subroutine



! ------------------------------------------------------------------------------
! SETDIFF routines


subroutine setdiff_real32 (arr1, arr2, n, diff, idx, assume_unique)
    !*  Compute the set difference A\B of sets A, B which are created
    !   from the unique elements of arr1 and arr2, respectively.

    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: arr1
        !*  Input array that defines A = set(arr1) in the operation A\B.
    real (PREC), intent(in), dimension(:) :: arr2
        !*  Input array that defined B = set(arr2) in the operation A\B.
    real (PREC), intent(out), dimension(:), optional :: diff
        !*  On exit, stores resulting set C = A\B that contains elements in arr1
        !   that are not in arr2. Elements are sorted in ascending order and
        !   duplicates are purged.
        !   If diff is not large enough not hold all elements of C, only the
        !   first size(diff) elements are stored and n > size(diff).

    real (PREC) :: val
#include "setdiff_impl.f90"
end subroutine


subroutine setdiff_real64 (arr1, arr2, n, diff, idx, assume_unique)
    !*  Compute the set difference A\B of sets A, B which are created
    !   from the unique elements of arr1 and arr2, respectively.

    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: arr1
        !*  Input array that defines A = set(arr1) in the operation A\B.
    real (PREC), intent(in), dimension(:) :: arr2
        !*  Input array that defined B = set(arr2) in the operation A\B.
    real (PREC), intent(out), dimension(:), optional :: diff
        !*  On exit, stores resulting set C = A\B that contains elements in arr1
        !   that are not in arr2. Elements are sorted in ascending order and
        !   duplicates are purged.
        !   If diff is not large enough not hold all elements of C, only the
        !   first size(diff) elements are stored and n > size(diff).

    real (PREC) :: val

#include "setdiff_impl.f90"
end subroutine


subroutine setdiff_int32 (arr1, arr2, n, diff, idx, assume_unique)
    !*  Compute the set difference A\B of sets A, B which are created
    !   from the unique elements of arr1 and arr2, respectively.

    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:) :: arr1
        !*  Input array that defines A = set(arr1) in the operation A\B.
    integer (INTSIZE), intent(in), dimension(:) :: arr2
        !*  Input array that defined B = set(arr2) in the operation A\B.
    integer (INTSIZE), intent(out), dimension(:), optional :: diff
        !*  On exit, stores resulting set C = A\B that contains elements in arr1
        !   that are not in arr2. Elements are sorted in ascending order and
        !   duplicates are purged.
        !   If diff is not large enough not hold all elements of C, only the
        !   first size(diff) elements are stored and n > size(diff).

    integer (INTSIZE) :: val

#include "setdiff_impl.f90"
end subroutine

!-------------------------------------------------------------------------------
! INTERSECT routines

subroutine intersect_real32 (arr1, arr2, res, n, assume_unique)
    !*  INTERSECT returns the unique, sorted elements that are present in
    !   both input arrays.

    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: arr1
        !*  First input array
    real (PREC), intent(in), dimension(:) :: arr2
        !*  Second input array
    real (PREC), intent(out), dimension(:) :: res
        !*  Intersection of input arrays, sorted in ascending order

    real (PREC), dimension(:), allocatable :: work

#include "intersect_impl.f90"
end subroutine


subroutine intersect_real64 (arr1, arr2, res, n, assume_unique)
    !*  INTERSECT returns the unique, sorted elements that are present in
    !   both input arrays.


    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: arr1
        !*  First input array
    real (PREC), intent(in), dimension(:) :: arr2
        !*  Second input array
    real (PREC), intent(out), dimension(:) :: res
        !*  Intersection of input arrays, sorted in ascending order

    real (PREC), dimension(:), allocatable :: work

#include "intersect_impl.f90"
end subroutine


subroutine intersect_int32 (arr1, arr2, res, n, assume_unique)
    !*  INTERSECT returns the unique, sorted elements that are present in
    !   both input arrays.


    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:) :: arr1
        !*  First input array
    integer (INTSIZE), intent(in), dimension(:) :: arr2
        !*  Second input array
    integer (INTSIZE), intent(out), dimension(:) :: res
        !*  Intersection of input arrays, sorted in ascending order

    integer (INTSIZE), dimension(:), allocatable :: work

#include "intersect_impl.f90"
end subroutine


!-------------------------------------------------------------------------------
! UNION routines


subroutine union_real32 (arr1, arr2, res, n, assume_unique)
    !*  INTERSECT returns the unique, sorted elements that are present in
    !   both input arrays.

    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: arr1
        !*  First input array
    real (PREC), intent(in), dimension(:) :: arr2
        !*  Second input array
    real (PREC), intent(out), dimension(:) :: res
        !*  Intersection of input arrays, sorted in ascending order

    real (PREC), dimension(:), allocatable :: work

#include "union_impl.f90"
end subroutine


subroutine union_real64 (arr1, arr2, res, n, assume_unique)
    !*  INTERSECT returns the unique, sorted elements that are present in
    !   both input arrays.


    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: arr1
        !*  First input array
    real (PREC), intent(in), dimension(:) :: arr2
        !*  Second input array
    real (PREC), intent(out), dimension(:) :: res
        !*  Intersection of input arrays, sorted in ascending order

    real (PREC), dimension(:), allocatable :: work

#include "union_impl.f90"
end subroutine


subroutine union_int32 (arr1, arr2, res, n, assume_unique)
    !*  INTERSECT returns the unique, sorted elements that are present in
    !   both input arrays.


    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:) :: arr1
        !*  First input array
    integer (INTSIZE), intent(in), dimension(:) :: arr2
        !*  Second input array
    integer (INTSIZE), intent(out), dimension(:) :: res
        !*  Intersection of input arrays, sorted in ascending order

    integer (INTSIZE), dimension(:), allocatable :: work

#include "union_impl.f90"
end subroutine



!-------------------------------------------------------------------------------
! SET_IN routines

pure function set_in_int32 (needle, haystack) result(res)
    !*  SET_IN operator returns .TRUE. whenever a given element
    !   is in a given set.
    !   No attempt has been made to optimize this routine, it is assumed
    !   that it's called only for small set arrays.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: needle
    integer (INTSIZE), intent(in), dimension(:) :: haystack
    logical :: res

    integer :: i

    res = .false.
    do i = 1, size(haystack)
        if (needle == haystack(i)) then
            res = .true.
            return
        end if
    end do

end function


pure function set_in_int64 (needle, haystack) result(res)
    !*  SET_IN operator returns .TRUE. whenever a given element
    !   is in a given set.
    !   No attempt has been made to optimize this routine, it is assumed
    !   that it's called only for small set arrays.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in) :: needle
    integer (INTSIZE), intent(in), dimension(:) :: haystack
    logical :: res

    integer :: i

    res = .false.
    do i = 1, size(haystack)
        if (needle == haystack(i)) then
            res = .true.
            return
        end if
    end do

end function



pure subroutine split_uniform_int32 (n, nchunks, ifrom, chunk_size, status)
    !*  SPLIT_UNIFORM partitions a sequence of size N into a given number
    !   of chunks which are of approximately uniform size.
    integer, intent(in) :: n
        !*  Number of elements to partition
    integer, intent(in) :: nchunks
        !*  Number of chunks (partitions) to create
    integer, intent(out), dimension(:) :: ifrom
        !*  Array to store the starting index of each chunk
    integer, intent(out), dimension(:) :: chunk_size
        !*  Array to store the number of elements in each chunk
    type (status_t), intent(out), optional :: status
        !*  Optional exit code

    integer :: k, i, ioffset
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'SPLIT_UNIFORM'

    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    call check_cond (n > 0, NAME, 'N: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (nchunks > 0, NAME, 'NCHUNKS: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (size(ifrom) /= size(chunk_size) .or. (size(ifrom) /= nchunks)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! --- Implementation ---

    if (n >= nchunks) then
        chunk_size = n / nchunks
        k = modulo (n, nchunks)
        chunk_size(1:k) = chunk_size(1:k) + 1
    else
        chunk_size = 0
        chunk_size(1:n) = 1
    end if

    ioffset = 1
    do i = 1, nchunks
        ifrom(i) = ioffset
        ioffset = ioffset + chunk_size(i)
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine



end module
