program test_nf_arrays_setops
    !*   Unit tests for array set routines.

    use, intrinsic :: iso_fortran_env
    use numfort_arrays

    use corelib_testing, only: test_suite, test_case
    use corelib_common, only: str

    implicit none

    integer, parameter :: PREC = real64

    call test_all()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_arrays_setops unit tests")

    call test_setdiff (tests)
    call test_intersect (tests)

    ! print test statistics
    call tests%print ()

end subroutine



subroutine test_setdiff (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: set_a, set_b, diff
    integer, dimension(:), allocatable :: idx
    integer :: n

    tc => tests%add_test ("setdiff routine")

    ! 1. Degenerate inputs A, B
    allocate (set_a(0), set_b(0), diff(0), idx(0))
    call setdiff (set_a, set_b, n, diff, idx)
    call tc%assert_true (n == 0, "Zero-size input sets A, B")
    deallocate (set_a, set_b, diff, idx)

    ! 2. Degenerate input B, some set A
    allocate (set_a(1), set_b(0), diff(1), idx(1))
    set_a = 1.0
    call setdiff (set_a, set_b, n, diff, idx)
    call tc%assert_true (n == 1 .and. all(set_a == diff) .and. all(idx == 1), &
        "Singleton A, zero-size B")
    deallocate (set_a, set_b, diff, idx)

    ! Degenerate set B, non-sorted set A with duplicates
    allocate (set_a(5), set_b(0), diff(5), idx(5))
    set_a(:) = [real (PREC) :: 1, 1, 5, 4, 7]
    idx = 0
    diff = 0
    call setdiff (set_a, set_b, n, diff, idx)
    ! unclear whether the underlying routine in ORDERPACK always
    ! returns index of last occurence of duplicate entries
    call tc%assert_true (n == 4 .and. all(diff(1:4) == [1,4,5,7]) .and. &
        all(idx(2:4) == [4,3,5]) .and. (idx(1) == 1 .or. idx(1) == 2), &
        "Zero-size B, unsorted A with duplicates")
    deallocate (set_a, set_b, diff, idx)

    ! 3. Degenerate input A, some set B
    allocate (set_a(0), set_b(1), diff(1), idx(1))
    set_b = 1.0
    call setdiff (set_a, set_b, n, diff, idx)
    call tc%assert_true (n == 0, "Zero-size A, singleton B")
    deallocate (set_a, set_b, diff, idx)


    ! 4. Identical non-degenerate set A, B
    allocate (set_a(5), set_b(5), diff(5), idx(5))
    call linspace (set_a, 1.0_PREC, 5.0_PREC)
    set_b(:) = set_a
    call setdiff (set_a, set_b, n, diff, idx)
    call tc%assert_true (n == 0, "Non-degenerate identical A, B")
    deallocate (set_a, set_b, diff, idx)


    ! 5. B is strict superset of A
    allocate (set_a(3), set_b(5), diff(5), idx(5))
    call linspace (set_b, 1.0_PREC, 5.0_PREC)
    set_a(:) = set_b(2:4)
    call setdiff (set_a, set_b, n, diff, idx)
    call tc%assert_true (n == 0, "Non-degenerate A, B; A strict subset of B")
    deallocate (set_a, set_b, diff, idx)


    ! 6. Some example that yields non-empty A\B
    allocate (set_a(5), set_b(3), diff(5), idx(5))
    diff = 0
    idx = 0
    call linspace (set_a, 1.0_PREC, 5.0_PREC)
    ! Create non-sorted input array
    set_a([1, 5]) = set_a([5, 1])
    set_b(:) = set_a(2:4)
    call setdiff (set_a, set_b, n, diff, idx)
    call tc%assert_true (n == 2 .and. all(diff(1:2)==[1,5]) .and. &
        all(idx(1:2) == [5,1]), &
        "Non-generate A, B; B strict subset of A")

    ! 7. Test with assume_unique = .true. using the same input arrays
    diff = 0
    idx = 0
    n = 0
    call setdiff (set_a, set_b, n, diff, idx, assume_unique=.true.)
    call tc%assert_true (n == 2 .and. all(diff(1:2)==[1,5]) .and. &
        all(idx(1:2) == [5,1]), &
        "Non-degenerate A, B; B strict subset of A; assumed unique")
    deallocate (set_a, set_b, diff, idx)


    ! 8. Test with B strict subset of A, B contains initial sorted elements
    ! of A.
    allocate (set_a(5), set_b(3), diff(2), idx(2))
    diff = 0
    idx = 0
    call linspace (set_a, 1.0_PREC, 5.0_PREC)
    set_b(1:3) = set_a(1:3)
    call setdiff (set_a, set_b, n, diff, idx, assume_unique=.true.)
    call tc%assert_true (n == 2 .and. all(diff(1:2)==[4,5]) .and. &
        all(idx(1:2)==[4,5]), &
        "B strict subset of A; B contains initial elems of A")

    ! Test with B containing last elements of A
    diff = 0
    idx = 0
    n = 0
    set_b(1:3) = set_a(3:5)
    call setdiff (set_a, set_b, n, diff, idx, assume_unique=.true.)
    call tc%assert_true (n == 2 .and. all(diff(1:2)==[1,2]) .and. &
        all(idx(1:2)==[1,2]), &
        "B strict subset of A; B contains last elemts of A")



end subroutine



subroutine test_intersect (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: set_a, set_b, res
    integer :: n, i

    tc => tests%add_test ("set INTERSECT routine: C = A intersect B")


    ! 1. Both operands are empty sets
    allocate (set_a(0), set_b(0), res(0))
    n = -1
    call intersect (set_a, set_b, res, n)
    call tc%assert_true (n == 0, "A, B empty sets")
    deallocate (set_a, set_b, res)

    ! 2. First operand empty, second operand non-empty
    allocate (set_a(0), set_b(10), res(10))
    set_b = [(i, i = 1, 10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    call tc%assert_true (n == 0, "A empty, B non-empty")
    deallocate (set_a, set_b, res)

    ! 3. Second operand empty, first non-empty
    allocate (set_a(10), set_b(0), res(0))
    set_a = [(i, i = 1, 10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    call tc%assert_true (n == 0, "A non-empty, B empty")
    deallocate (set_a, set_b, res)

    ! 4. Nonempty operands, empty intersection
    allocate (set_a(5), set_b(5), res(0))
    set_a = [(i, i=1,5)]
    set_b = [(i, i=6,10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    call tc%assert_true (n == 0, "A, B non-empty, empty intersection")
    deallocate (set_a, set_b, res)

    ! 5. A strict subset of B
    allocate (set_a(5), set_b(10), res(5))
    set_a = [(i, i = 1, 5)]
    set_b = [(i, i = 1, 10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    call tc%assert_true (n == 5 .and. all(res(1:n) == set_a), &
        "A, B nonempty, C == A")
    deallocate (set_a, set_b, res)

    ! 6. B strict subset of A
    allocate (set_a(10), set_b(5), res(5))
    set_a = [(i, i = 1, 10)]
    set_b = [(i, i = 6, 10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    call tc%assert_true (n == 5 .and. all(res(1:n) == set_b), &
        "A, B nonempty, C == B")
    deallocate (set_a, set_b, res)

    ! 7. Test some random non-empty intersection
    allocate (set_a(5), set_b(5), res(3))
    set_a = [(i, i = 1, 5)]
    set_b = [2, 3, 5, 7, 10]
    n = -1
    call intersect (set_a, set_b, res, n)
    call tc%assert_true (n == 3 .and. all(res(1:n) == [2,3,5]), &
        "Non-empty A, B; non-empty intersection")
    deallocate (set_a, set_b, res)

end subroutine



end program
