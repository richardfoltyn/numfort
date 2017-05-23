module numfort_arrays_setops

    use, intrinsic :: iso_fortran_env
    use m_unirnk, only: unirnk
    use numfort_arrays_sort

    implicit none

    private

    public :: unique, setdiff

    interface unique
        module procedure unique_real32, unique_real64, unique_int32
    end interface

    interface setdiff
        module procedure setdiff_real32, setdiff_real64, setdiff_int32
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

end module