program test_arrays_vander
    !*  Unit tests for Vandermonde matrix creation routines.

    use, intrinsic :: iso_fortran_env
    use numfort_arrays, only: vander
    use numfort_common_testing, only: all_close

    use fcore_testing, only: test_suite, test_case
    use fcore_common, only: str

    implicit none

    integer, parameter :: PREC = real64

    call test_all()

    contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("Vandermonde matrix creation unit tests")

    call test_scalar (tests)
    call test_array (tests)

    call tests%print ()

end subroutine


subroutine test_scalar (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: xp, xp0
    real (PREC) :: x
    integer :: np, j

    tc => tests%add_test ("Vandermonde from scalar argument")

    np = 5
    allocate (xp(np), xp0(np))

    x = 0.0
    call vander (x, xp)
    call tc%assert_true (all(xp(:np-1)==0.0) .and. xp(np) == 1.0, "Scalar input, x=0.0")

    x = 1.0
    call vander (x, xp)
    call tc%assert_true (all(xp==1.0), "Scalar input, x=1.0")

    x = 0.123
    do j = 0, np-1
        xp0(j+1) = x ** (np-1-j)
    end do
    call vander (x, xp)
    call tc%assert_true (all_close(xp, xp0), "Scalar input, default order")

    x = 12.345
    do j = 0, np-1
        xp0(j+1) = x ** j
    end do
    call vander (x, xp, increasing=.true.)
    call tc%assert_true (all_close(xp, xp0), "Scalar input, increasing order")

end subroutine


subroutine test_array (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: x
    real (PREC), dimension(:,:), allocatable :: xp, xp0

    integer :: np, nx, j

    tc => tests%add_test ("Vandermonde from vector argument")

    ! Test with degenerate input arrays. Just verifies that routines don't
    ! crash
    allocate (x(0), xp(0, 0))
    call vander (x, xp)
    deallocate (x, xp)

    allocate (x(1), xp(0, 1))
    call vander (x, xp)
    deallocate (x, xp)

    allocate (x(1), xp(1, 0))
    call vander (x, xp)
    deallocate (x, xp)

    np = 5
    allocate (x(1), xp(2, np), xp0(1, np))
    x(1) = 2.0
    do j = 0, np-1
        xp0(:,j+1) = x**(np-1-j)
    end do
    call vander (x, xp)
    call tc%assert_true (all_close (xp(1,:), xp(1,:)), &
        "Vector input, size(x) < size(xp,1)")
    deallocate (x, xp,  xp0)

    np = 5
    allocate (x(2), xp(1,np), xp0(1,np))
    x(1) = 123.456
    do j = 0, np-1
        xp0(:,j+1) = x(1)**(np-1-j)
    end do
    call vander (x, xp)
    call tc%assert_true (all_close (xp(1,:), xp(1,:)), &
        "Vector input, size(x) > size(xp,1)")
    deallocate (x, xp,  xp0)

    ! Test with some non-degenerate input vector
    nx = 3
    np = 4

    ! 1. exponents in decreasing order (default)
    allocate (x(nx), xp(nx,np), xp0(nx,np))
    x(:) = [0.1d0, 1.5d0, 10.123d0]
    do j = 0, np-1
        xp0(:,j+1) = x**(np-1-j)
    end do
    call vander (x, xp)
    call tc%assert_true (all(xp0==xp), "Vector input, default order")

    ! 2. Exponents in increasing order
    do j = 0, np-1
        xp0(:,j+1) = x**j
    end do
    call vander (x, xp, increasing=.true.)
    call tc%assert_true (all(xp0==xp), "Vector input, increasing order")

    deallocate (x, xp)


end subroutine

end program
