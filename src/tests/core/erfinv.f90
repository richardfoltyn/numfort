

program numfort_test_core_erfinv
    !*  Unit tests for inverse ERF function (ERFINV)

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays, only: linspace
    use numfort_core, only: erfinv
    use numfort_common_testing, only: all_close


    use fcore_testing
    use fcore_strings

    implicit none


    integer, parameter :: PREC = real64


    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("ERFINV unit tests")

    call test_erfinv (tests)

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

end program
