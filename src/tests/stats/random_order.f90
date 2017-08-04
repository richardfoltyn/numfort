program test_numfort_stats_random_order

    use iso_fortran_env

    use fcore_common
    use fcore_testing
    use numfort_stats, only: random_order

    implicit none

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_stats_combinatorics unit tests")

    call test_random_order (tests)
    ! call test_sub2ind (tests)

    ! print test statistics
    call tests%print ()

end subroutine


subroutine test_random_order (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, dimension(:), allocatable :: x
    integer :: n, low
    logical correct

    tc => tests%add_test ("random_order test cases")

    n = 100
    ! call with default low parameter
    allocate (x(n))
    call random_order (x)
    call random_order_verify (x, 1, correct)

    call tc%assert_true (correct, &
        "argument x(100), default low parameter")
    deallocate (x)

    ! Check with non-default positive argument
    allocate (x(n))
    low = 11
    call random_order (x, low=low)
    call random_order_verify (x, low, correct)
    call tc%assert_true (correct, "argument x(100), low=" // str(low))
    deallocate (x)

    ! Check with non-default negative argument
    allocate (x(n))
    low = -1713
    call random_order (x, low=low)
    call random_order_verify (x, low, correct)
    call tc%assert_true (correct, "argument x(100), low=" // str(low))
    deallocate (x)

    ! Check with degenerate size(x) = 1
    n = 1
    allocate (x(n))
    low = 1
    call random_order (x, low=low)
    call random_order_verify (x, low, correct)
    call tc%assert_true (correct, "argument x(1), low=" // str(low))
    deallocate (x)

    ! Check with degenerate size(x) = 0
    n = 0
    allocate (x(n))
    low = 1
    call random_order (x, low=low)
    ! nothing to assert, just check that it does not crash


end subroutine

subroutine random_order_verify (x, low, correct)
    integer, intent(in), dimension(:) :: x
    integer, intent(in) :: low
    logical, intent(out) :: correct

    integer, dimension(:), allocatable :: counter
    integer :: i, j, imin, imax

    allocate (counter(size(x)), source = 0)

    ! check that all integers {low,...,low+n-1} are present exactly once
    do i = 1, size(x)
        j = x(i) - low + 1
        counter(j) = counter(j) + 1
    end do

    imin = minval(x)
    imax = maxval(x)

    correct = all(counter == 1) .and. (imin == low) .and. (imax == low + size(x) - 1)
end subroutine


end program
