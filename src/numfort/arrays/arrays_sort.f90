module numfort_arrays_sort

    use, intrinsic :: iso_fortran_env
    ! mergesort implementation from ORDERPACK
    use m_mrgrnk, only: mrgrnk
    use numfort_common

    implicit none
    private

    !>  Return the indices that would sort an array, ie the rank of
    !   each array element.
    interface argsort
        module procedure argsort_real32, argsort_real64, argsort_int32
    end interface

    public :: argsort

contains

subroutine argsort_real32 (arr, rnk, algorithm, status)
    !*  Return the indices that would sort an array, ie the rank of
    !   each array element.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:) :: arr
        !!  Array to sort

    include "include/argsort_impl.f90"
end subroutine

subroutine argsort_real64 (arr, rnk, algorithm, status)
    !*  Return the indices that would sort an array, ie the rank of
    !   each array element.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:) :: arr
        !!  Array to sort

    include "include/argsort_impl.f90"
end subroutine

subroutine argsort_int32 (arr, rnk, algorithm, status)
    !*  Return the indices that would sort an array, ie the rank of
    !   each array element.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:) :: arr
        !!  Array to sort

    include "include/argsort_impl.f90"
end subroutine

end module
