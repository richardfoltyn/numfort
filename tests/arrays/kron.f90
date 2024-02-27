program test_arrays_kron
    !*  Unit tests for Vandermonde matrix creation routines.

    use, intrinsic :: iso_fortran_env
    use numfort_arrays, only: kron
    use numfort_common_testing, only: all_close

    use fcore_testing, only: test_suite, test_case
    use fcore_common, only: str

    implicit none

    integer, parameter :: PREC = real64

    call test_all()

    contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("Kronecker product unit tests")

    call test_kron (tests)

    call tests%print ()

end subroutine


subroutine test_kron (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(:,:), allocatable :: m1, m2, mm, mm_ok

    tc => tests%add_test ("Kronecker product with valid arguments")

    allocate (m1(2,2), m2(2,2))
    m1(:,:) = reshape([1,2,3,4], shape=[2,2])
    m2(:,:) = reshape([1,0,0,1], shape=[2,2])

    allocate (mm(4,4), mm_ok(4,4))
    call kron (m1, m2, mm)

    mm_ok(1:2,1:2) = m2
    mm_ok(3:4,1:2) = m1(2,1) * m2
    mm_ok(1:2,3:4) = m1(1,2) * m2
    mm_ok(3:4,3:4) = m1(2,2) * m2

    call tc%assert_true (all_close (mm, mm_ok, atol=0.0_PREC, rtol=0.0_PREC), &
        "Kronecker product of two 2x2 matrices")


end subroutine


end program
