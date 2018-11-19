

program test_polynomial_herme
    !*  Unit tests for "probabilist" Hermite polynomial routines.

    use, intrinsic :: iso_fortran_env

    use numfort_polynomial_hermite_e
    use numfort_common
    use numfort_common_testing
    use numfort_io

    use fcore_testing, only: test_suite, test_case
    use fcore_strings


    implicit none

    integer, parameter :: PREC = real64

    call test_all ()


    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("Hermite_e polynomial unit tests")

    call test_quadrature (tests)

    call tests%print ()

end subroutine


subroutine test_quadrature (tests)
    !*  TEST_QUADRATURE runs unit tests for routines computing sample
    !   points and weights for Gauss-Hermite quadrature.
    !
    !   Note that the code expects the CSV output generated by numpy's
    !   hermegauss() routine to be placed in the data/ subdirectory
    !   underneath the working directory from which the unit test
    !   executable is run.
    !
    !   Note: the input data is generated by running the Python script
    !   data/generate_quad_data.py
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (str) :: msg

    real (PREC), dimension(:), allocatable :: x, w
    integer :: i, n, k, ifrom, ito
    type (status_t) :: status
    integer, dimension(:), allocatable :: degrees, iwork
    real (PREC), dimension(:), allocatable :: x_np, w_np
    real (PREC), dimension(:,:), allocatable :: dat
    logical :: x_ok, w_ok

    tc => tests%add_test ("Gaussian quadrature")

    ! === Test with invalid inputs ===
    ! 1. zero-size input
    n = 0
    allocate (x(n), w(n))
    status = NF_STATUS_UNDEFINED
    call hermegauss (x, w, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, "Zero-size input arrays")
    deallocate (x, w)

    ! 2. test with different-size input arrays
    allocate (x(1), w(2))
    status = NF_STATUS_UNDEFINED
    call hermegauss (x, w, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Different-size input arrays")
    deallocate (x, w)

    ! === Compare to numpy output ===

    ! Read in list of number of quadrature points (=polynomial degrees)
    allocate (iwork(100), source=0)
    call read_fixed (path='data/hermegauss_n.csv', fmt='(*(i3,:,/))', dat=iwork)

    ! drop trailing zeros
    k = count (iwork > 0)
    allocate (degrees(k), source=iwork(1:k))

    allocate (dat(3, sum(degrees)), source=0.0_PREC)
    ! Read in actual quadrature data
    call read_fixed (path='data/hermegauss_data.csv', fmt='(*(f3.0, 2(f50.30)))', &
        dat=dat, transform='transpose', status=status)

    do i = 1, k
        n = degrees(i)
        allocate (x_np(n), w_np(n))
        allocate (x(n), w(n))

        status =  NF_STATUS_UNDEFINED
        call hermegauss (x, w, status=status)

        ifrom = sum(degrees(1:i-1)) + 1
        ito = sum(degrees(1:i))

        x_np(:) = dat(2, ifrom:ito)
        w_np(:) = dat(3, ifrom:ito)

        x_ok = all_close (x, x_np, rtol=0.0_PREC, atol=1.0e-12_PREC)
        w_ok = all_close (w, w_np, rtol=0.0_PREC, atol=1.0e-12_PREC)

        msg = 'Equals numpy solution for ' // str(n, 'i3') // ' quadrature nodes'
        call tc%assert_true (x_ok .and. w_ok .and. status == NF_STATUS_OK, &
            msg)

        deallocate (x, w, x_np, w_np)
    end do

end subroutine


end program