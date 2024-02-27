
program test_arrays_delete
    !*  Module contains unit test routines for the array DELETE routine.

    use, intrinsic :: iso_fortran_env

    use numfort_common_status
    use numfort_arrays, only: delete

    use fcore_testing, only : test_suite, test_case

    implicit none

    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("numfort_arrays DELETE unit tests")

    ! call test_inplace (tests)
    call test_1d_copy (tests)
    call test_2d_copy (tests)

    call tests%print ()

end subroutine


subroutine test_1d_copy (tests)
    !*  Unit tests for deleting elements and storing result in a copy.
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: PREC = real32
    integer, parameter :: na = 5, nr = 5
    real (PREC) :: arr(na), res(nr)
    type (status_t) :: status
    integer :: i, n

    tc => tests%add_test ("DELETE elements from 1d-array, store result in copy")
    arr = [(i, i=1,na)]

    ! Test invalid indices
    res = -1
    status = NF_STATUS_UNDEFINED
    call delete (arr, 0, res, n=n, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG .and. all(res==-1) &
        .and. n == 0, &
        "DELETE called with IDX < 1")

    status = NF_STATUS_UNDEFINED
    call delete (arr, na+1, res, n=n, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG .and. all(res==-1) &
        .and. n == 0, &
        "DELETE called with IDX > size(ARR)")

    ! Test with empty index array
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [integer ::], res, n=n, status=status)
    call tc%assert_true (all(arr==res) .and. status == NF_STATUS_OK .and. n==na, &
        "DELETE called with empty IDX array")

    ! Test deleting element at beginning of arr
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, 1, res, n=n, status=status)
    call tc%assert_true (all(arr(2:)==res(1:na-1)) .and. status == NF_STATUS_OK, &
        "DELETE called with IDX=1")

    ! Test deleting element at end of ARR
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, na, res, n=n, status=status)
    call tc%assert_true (all(arr(1:na-1)==res(1:na-1)) &
        .and. status == NF_STATUS_OK .and. n==(na-1), &
        "DELETE called with IDX=size(ARR)")

    ! Test deleting multile interior elements
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [2,3], res, n=n, status=status)
    call tc%assert_true (arr(1)==res(1) .and. all(arr(4:5)==res(2:3)) &
        .and. status == NF_STATUS_OK .and. n==3, &
        "DELETE called with IDX=[2,3]")

    ! Test with RES array that is not large enough
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [2,3], res(1:2), n=n, status=status)
    call tc%assert_true (arr(1)==res(1) .and. arr(4)==res(2) .and. &
        all(res(3:)==-1) .and. status == NF_STATUS_OK .and. n==2, &
        "DELETE called with IDX=[2,3] and OUT not large enough")

    ! Test with unsorted index array
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2], res, n=n, status=status)
    call tc%assert_true (arr(1)==res(1) .and. all(arr(4:5)==res(2:3)) &
        .and. status == NF_STATUS_OK .and. n==3, &
        "DELETE called with IDX=[3,2] (unsorted)")

    ! Test with non-unique index array
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2,2], res, n=n, status=status)
    call tc%assert_true (arr(1)==res(1) .and. all(arr(4:5)==res(2:3)) &
        .and. status == NF_STATUS_OK .and. n==3, &
        "DELETE called with IDX=[3,2,2] (unsorted & non-unique)")

    ! Test without passing any optional parameters
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2,2], res)
    call tc%assert_true (arr(1)==res(1) .and. all(arr(4:5)==res(2:3)), &
        "DELETE called with IDX=[3,2,2], no optional args")

    ! Test with non-contiguous block of indices to be deleted
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2,5], res)
    call tc%assert_true (arr(1)==res(1) .and. arr(4)==res(2), &
        "DELETE called with IDX=[3,2,5], no optional args")
    
    ! Test with IDX such that all elements are deleted
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [(i,i=1,na)], res, n=n)
    call tc%assert_true (all(res==-1) .and. n==0, &
        "DELETE called with IDX=[ALL]")

end subroutine


subroutine test_2d_copy (tests)
    !*  Unit tests for deleting elements and storing result in a copy.
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: PREC = real32
    integer, parameter :: na = 5, nr = 5
    real (PREC) :: arr(na,na), res(nr,nr)
    type (status_t) :: status
    integer :: i, n

    tc => tests%add_test ("DELETE elements from 2d-array, store result in copy")
    arr = reshape([(i, i=1,na*na)], shape=[na,na])

    ! Test invalid indices
    res = -1
    status = NF_STATUS_UNDEFINED
    call delete (arr, 0, 1, res, n=n, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG .and. all(res==-1) &
        .and. n == 0, &
        "DELETE called with IDX < 1")

    status = NF_STATUS_UNDEFINED
    call delete (arr, na+1, 1, res, n=n, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG .and. all(res==-1) &
        .and. n == 0, &
        "DELETE called with IDX > size(ARR)")

    ! Test with empty index array
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [integer ::], 1, res, n=n, status=status)
    call tc%assert_true (all(arr==res) .and. status == NF_STATUS_OK .and. n==na, &
        "DELETE called with empty IDX array")

    ! Test with result array that has too few rows
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [1], 1, res(1:3,:), n=n, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG .and. all (res==-1) &
        .and. n == 0, &
        "DELETE called with OUT that has too few rows")

    ! Test with result array that has too few columns
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [1], 1, res(1:4,1:4), n=n, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG .and. all (res==-1) &
        .and. n == 0, &
        "DELETE called with OUT that has too few columns")

    ! Test deleting element at beginning of arr
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, 1, 1, res, n=n, status=status)
    call tc%assert_true (all(arr(2:,:)==res(1:na-1,:)) &
        .and. status == NF_STATUS_OK .and. n==(na-1), &
        "DELETE called with IDX=1, DIM=1")

    ! Test deleting element at end of ARR
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, na, 1, res, n=n, status=status)
    call tc%assert_true (all(arr(1:na-1,:)==res(1:na-1,:)) &
        .and. status == NF_STATUS_OK .and. n==(na-1), &
        "DELETE called with IDX=size(ARR,1), DIM=1")

    ! Test deleting element at beginning of arr
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, 1, 2, res, n=n, status=status)
    call tc%assert_true (all(arr(:,2:)==res(:,1:na-1)) &
        .and. status == NF_STATUS_OK .and. n==(na-1), &
        "DELETE called with IDX=1, DIM=2")

    ! Test deleting element at end of ARR
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, na, 2, res, n=n, status=status)
    call tc%assert_true (all(arr(:,1:na-1)==res(:,1:na-1)) &
        .and. status == NF_STATUS_OK .and. n==(na-1), &
        "DELETE called with IDX=size(ARR,1), DIM=2")

    ! Test deleting multile interior elements
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [2,3], 1, res, n=n, status=status)
    call tc%assert_true (all(arr(1,:)==res(1,:)) .and. all(arr(4:5,:)==res(2:3,:)) &
        .and. status == NF_STATUS_OK .and. n==3, &
        "DELETE called with IDX=[2,3], DIM=1")

    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [2,3], 2, res, n=n, status=status)
    call tc%assert_true (all(arr(:,1)==res(:,1)) .and. all(arr(:,4:5)==res(:,2:3)) &
        .and. status == NF_STATUS_OK .and. n==3, &
        "DELETE called with IDX=[2,3], DIM=2")

    ! Test with unsorted index array
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2], 1, res, n=n, status=status)
    call tc%assert_true (all(arr(1,:)==res(1,:)) .and. all(arr(4:5,:)==res(2:3,:)) &
        .and. status == NF_STATUS_OK .and. n==3, &
        "DELETE called with IDX=[3,2] (unsorted), DIM=1")

    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2], 2, res, n=n, status=status)
    call tc%assert_true (all(arr(:,1)==res(:,1)) .and. all(arr(:,4:5)==res(:,2:3)) &
        .and. status == NF_STATUS_OK .and. n==3, &
        "DELETE called with IDX=[3,2] (unsorted), DIM=2")

    ! Test with non-unique index array
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2,2], 1, res, n=n, status=status)
    call tc%assert_true (all(arr(1,:)==res(1,:)) .and. all(arr(4:5,:)==res(2:3,:)) &
        .and. status == NF_STATUS_OK .and. n==3, &
        "DELETE called with IDX=[3,2,2] (unsorted & non-unique), DIM=1")

    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2,2], 2, res, n=n, status=status)
    call tc%assert_true (all(arr(:,1)==res(:,1)) .and. all(arr(:,4:5)==res(:,2:3)) &
        .and. status == NF_STATUS_OK .and. n==3, &
        "DELETE called with IDX=[3,2,2] (unsorted & non-unique), DIM=2")

    ! Test without passing any optional parameters
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2,2], 1, res)
    call tc%assert_true (all(arr(1,:)==res(1,:)) .and. all(arr(4:5,:)==res(2:3,:)), &
        "DELETE called with IDX=[3,2,2], DIM=1, no optional args")

    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2,2], 2, res)
    call tc%assert_true (all(arr(:,1)==res(:,1)) .and. all(arr(:,4:5)==res(:,2:3)), &
        "DELETE called with IDX=[3,2,2], DIM=2, no optional args")

    ! Test with non-contiguous block of indices to be deleted
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2,5], 1, res)
    call tc%assert_true (all(arr(1,:)==res(1,:)) .and. all(arr(4,:)==res(2,:)), &
        "DELETE called with IDX=[3,2,5], DIM=1, no optional args")

    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [3,2,5], 2, res)
    call tc%assert_true (all(arr(:,1)==res(:,1)) .and. all(arr(:,4)==res(:,2)), &
        "DELETE called with IDX=[3,2,5], DIM=2, no optional args")
    
    ! Test with IDX such that all elements are deleted
    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [(i,i=1,na)], 1, res, n=n)
    call tc%assert_true (all(res==-1) .and. n==0, &
        "DELETE called with IDX=[ALL], DIM=1")

    status = NF_STATUS_UNDEFINED
    res = -1
    call delete (arr, [(i,i=1,na)], 2, res, n=n)
    call tc%assert_true (all(res==-1) .and. n==0, &
        "DELETE called with IDX=[ALL], DIM=2")

end subroutine


end program
