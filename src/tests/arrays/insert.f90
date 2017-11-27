
program test_arrays_insert
    !*  Module contains unit test routines for the array INSERT routine.

    use, intrinsic :: iso_fortran_env
    use numfort_arrays, only: insert

    use fcore_testing

    implicit none

    integer, parameter :: PREC = real64


    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("numfort_arrays INSERT unit tests")

    call test_inplace (tests)

    call tests%print ()

end subroutine


subroutine test_inplace (tests)
    !*  Unit tests for in-place insertion into given argument array.

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: na = 5, nv = 3
    integer :: arr(na), val(nv), arr0(na), res(na)
    integer :: i

    tc => tests%add_test ("INSERT elements in place")

    arr = [(i, i=1,size(arr))]
    ! store copy of original array
    arr0 = arr
    val = [(10+i, i=1,size(val))]

    ! Invalid index
    call insert (arr, 0, val)
    call tc%assert_true (all(arr==arr0), "INSERT with idx=0")

    call insert (arr, size(arr)+2, val)
    call tc%assert_true (all(arr==arr0), "INSERT with idx >= size(arr)+1")

    ! Insert empty VAL array
    call insert (arr, 1, val(:0))
    call tc%assert_true (all(arr==arr0), "INSERT with empty VAL argument")

    ! Insert VAL at beginning
    call insert (arr, 1, val)
    call tc%assert_true (all(arr(1:nv)==val) .and. all(arr(4:)==arr0(:2)), &
        "Insert non-empty VAL at the beginning")

    ! Insert VAL at end
    arr = arr0
    call insert (arr, na-nv+1, val)
    call tc%assert_true (all(arr(:na-nv)==arr0(:na-nv)) .and. &
        all(arr(na-nv+1:)==val), &
        "Insert non-empty VAL at the end")

    ! Insert VAL in the interior
    arr = arr0
    call insert (arr, 2, val)
    call tc%assert_true (all(arr(2:2+nv-1)==val) .and. arr(1)==arr0(1) &
        .and. arr(na)==arr0(2), &
        "Insert non-empty VAL in the interior of ARR")

    ! Insert VAL such that it won't fit completely
    arr = arr0
    call insert (arr, 4, val)
    call tc%assert_true (all(arr(:3)==arr0(:3)) .and. all(arr(4:)==val(:2)), &
        "Insert VAL at end that won't fit completely")

    arr = arr0
    call insert (arr(1:2), 1, val)
    call tc%assert_true (all(arr(:2)==val(:2)) .and. all(arr(3:)==arr0(3:)), &
        "Insert VAL at the beginning that won't fit completely")

end subroutine

end program
