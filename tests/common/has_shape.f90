

program test_common_has_shape
    !*  Module contains unit test routines for the array INSERT routine.

    use, intrinsic :: iso_fortran_env
    use numfort_common, only: has_shape

    use fcore_testing
    use fcore_common, only: str, str_array, operator(//)

    implicit none

    integer, parameter :: PREC = real64


    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("numfort_common HAS_SHAPE unit tests")

    call test_1d (tests)
    call test_2d (tests)

    call tests%print ()

end subroutine


subroutine test_1d (tests)
    !*  Unit tests for in-place insertion into given argument array.

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: n = 2
    integer, dimension(n) :: arr
    logical :: res, res0
    integer :: i, j
    type (str) :: msg

    tc => tests%add_test ("1d-array tests")

    do i = 0, n+1
        do j = 0, n
            ! Test with scalar interface
            res = has_shape (arr(:j), i)
            res0 = size(arr(:j)) == i
            msg = "ARR of size " // str(j, "i0") // "; SHP=" // str(i, "i0")
            call tc%assert_true (res .eqv. res0, msg)

            ! Test with empty arrays, array interface
            res = has_shape (arr(:j), [i])
            res0 = all(shape(arr(:j)) == [i])
            msg = "ARR of size " // str(j, "i0") // "; SHP=[" // str(i, "i0") // "]"
            call tc%assert_true (res .eqv. res0, msg)
        end do
    end do

end subroutine


subroutine test_2d (tests)
    !*  Unit tests for in-place insertion into given argument array.

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: n = 2
    integer, dimension(n,n) :: arr
    integer, dimension(:), allocatable :: shp
    logical :: res, res0
    integer :: i, j, k, l
    type (str) :: msg

    tc => tests%add_test ("2d-array tests")

    do i = 0, n+1
        do j = 0, n+1
            do k = 0, n
                do l = 0, n
                    res = has_shape (arr(:k,:l), [i,j])
                    res0 = all(shape(arr(:k,:l)) == [i,j])
                    msg = "ARR of shape " // str_array(shape(arr(:k,:l)), "i0") &
                        // "; SHP=" // str_array([i,j], "i0")
                    call tc%assert_true (res .eqv. res0, msg)
                end do
            end do
        end do
    end do

    ! Check with incompatible shape parameter, should return .false.
    ! instead of crashing the application
    do i = 0, 5
        if (i == size(shape(arr))) cycle

        allocate (shp(i), source=0)
        res = has_shape (arr, shp)
        msg = '2d-array, shape argument of length ' // str(i)
        call tc%assert_false (res, msg)

        deallocate (shp)
    end do

end subroutine


end program
