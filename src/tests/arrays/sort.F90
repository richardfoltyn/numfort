


program test_numfort_arrays_sort

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common_status

    use numfort_common_testing

    use numfort_stats

    use fcore_strings
    use fcore_testing

    implicit none


    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("numfort_arrays sorting unit tests")

    call test_argsort_real32_int32 (tests)
    call test_argsort_real32_int64 (tests)
    call test_argsort_real64_int32 (tests)
    call test_argsort_real64_int64 (tests)

    call test_argsort_int32 (tests)
    call test_argsort_int64 (tests)

    call tests%print ()

end subroutine



subroutine test_argsort_real32_int32 (tests)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    tc => tests%add_test ('Unit test for ARGSORT (real32/int32')

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine



subroutine test_argsort_real32_int64 (tests)
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    tc => tests%add_test ('Unit test for ARGSORT (real32/int64')

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine



subroutine test_argsort_real64_int32 (tests)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    tc => tests%add_test ('Unit test for ARGSORT (real64/int32')

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine



subroutine test_argsort_real64_int64 (tests)
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int64

    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    tc => tests%add_test ('Unit test for ARGSORT (real64/int64')

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine



subroutine test_argsort_int32 (tests)
    integer, parameter :: INTSIZE = int32

    class (test_suite) :: tests

    class (test_case), pointer :: tc

    integer (INTSIZE), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    tc => tests%add_test ('Unit test for ARGSORT (int32')

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine



subroutine test_argsort_int64 (tests)
    integer, parameter :: INTSIZE = int64

    class (test_suite) :: tests

    class (test_case), pointer :: tc

    integer (INTSIZE), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    tc => tests%add_test ('Unit test for ARGSORT (int64')

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine

end
