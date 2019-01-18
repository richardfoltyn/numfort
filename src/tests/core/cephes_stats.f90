

program numfort_test_core_cephes_stats
    !*  Unit tests for selected stats-related functions ported from the
    !   the Cephes math library.

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays, only: linspace
    use numfort_core, only: erfinv, ndtr, ndtri
    use numfort_common_testing, only: all_close


    use fcore_testing
    use fcore_strings

    implicit none


    integer, parameter :: PREC = real64


    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("Unit tests for stats-related functions")

    call test_erfinv (tests)
    call test_ndtr (tests)

    call tests%print ()

end subroutine



subroutine test_erfinv (tests)
    !*  Test ERFINV inverse function by checking for approximate identity
    !   of ERFINV(ERF(x)) where ERF() is the built-in ERF function.
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC) :: xi, yi, x1, diff
    real (PREC), dimension(:), allocatable :: xx
    integer :: n, i
    real (PREC) :: tol = 1.0e-10_PREC
    type (str) :: msg

    tc => tests%add_test ("ERFINV unit tests")

    n = 20
    allocate (xx(20))
    call linspace (xx, 0.0_PREC, 3.5_PREC)

    do i = 1, size(xx)
        xi = xx(i)
        yi = erf (xi)
        x1 = erfinv (yi)
        diff = abs(xi-x1)
        msg = 'ERFINV(y) for y =' // str(yi, 'es10.4e2') &
            // '; abs(x0-x1) = ' // str(diff, 'es8.2e2')
        call tc%assert_true (abs(x1-xx(i)) < tol, msg)
    end do

end subroutine



subroutine test_ndtr (tests)
    !*  Unit tests for Gaussian CDF and inverse CDF routines NDTR and NDTRI.
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC) :: xi, yi, x1, y1, diff
    real (PREC), dimension(:), allocatable :: xx, yy
    integer :: n, i

    real (PREC) :: tol = 1.0e-10_PREC
    type (str) :: msg

    tc => tests%add_test ("NDTR/NDTRI unit tests")

    n = 20
    allocate (xx(20))
    call linspace (xx, -4.0_PREC, 4.0_PREC)

    do i = 1, size(xx)
        xi = xx(i)
        yi = ndtr (xi)
        x1 = ndtri (yi)
        diff = abs(x1-xi)
        msg = 'NDTRI(NDTR(x)) for x = ' // str(xi, 'es10.3e2') &
            // '; abs(x0-x1) = ' // str(diff, 'es8.2e2')
        call tc%assert_true (diff < tol, msg)
    end do

    allocate (yy(20))
    call linspace (yy, 1.0e-10_PREC, 1.0_PREC-1.0e-10_PREC)

    do i = 1, size(yy)
        yi = yy(i)
        xi = ndtri (yi)
        y1 = ndtr (xi)
        diff = abs(y1-yi)
        msg = 'NDTR(NDTRI(y)) for y = ' // str(yi, 'es10.3e2') &
            // '; abs(y0-y1) = ' // str(diff, 'es8.2e2')
        call tc%assert_true (diff < tol, msg)
    end do

    ! Test at some well-known values
    yi = 0.0
    xi = ndtri (yi)
    call tc%assert_true (ieee_class (xi) == IEEE_NEGATIVE_INF, &
        'NDTRI(0.0) == -inf')

    yi = 1.0
    xi = ndtri (yi)
    call tc%assert_true (ieee_class (xi) == IEEE_POSITIVE_INF, &
        'NDTRI(1.0) == inf')

    xi = ieee_value (xi, IEEE_NEGATIVE_INF)
    yi = ndtr (xi)
    call tc%assert_true (yi == 0.0_PREC, 'NDTR(-inf) == 0.0')

    xi = ieee_value (xi, IEEE_POSITIVE_INF)
    yi = ndtr (xi)
    call tc%assert_true (yi == 1.0_PREC, 'NDTR(inf) == 1.0')

    xi = 0.0
    yi = ndtr (xi)
    call tc%assert_true (abs(yi-0.5_PREC) < 1.0e-15_PREC, 'NDTR(0.0) == 0.5')

    yi = 0.5
    xi = ndtri (yi)
    call tc%assert_true (abs(xi) < 1.0e-15_PREC, 'NDTRI(0.5) == 0.0')

end subroutine


end program
