



program test_numfort_core_libc

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: linspace
    use numfort_core
    use numfort_common_testing

    use fcore_strings
    use fcore_testing

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ('Unit tests for special LOG/EXP functions')

    call test_logaddexp (tests)

    call tests%print ()

end subroutine


subroutine test_logaddexp (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: xx, yy, fxx, fxx1
    integer :: n

    tc => tests%add_test ('Unit tests for LOGADDEXP')

    ! Test for non-critical values, ie not for large negative values
    n = 100
    allocate (xx(n), yy(n), fxx(n), fxx1(n))

    call linspace(xx, -20.0_PREC, 10.0_PREC)
    call linspace(yy, -10.0_PREC, 15.0_PREC)

    fxx(:) = log(exp(xx) + exp(yy))
    fxx1(:) = logaddexp(xx, yy)

    call tc%assert_true (all_close (fxx1, fxx, atol=1.0e-10_PREC), &
        'Test LOGADDEXP(x,y) = LOG(EXP(X) + EXP(Y)) for non-critical values')

    fxx(:) = log(exp(yy) + exp(xx))
    fxx1(:) = logaddexp(yy, xx)
    call tc%assert_true (all_close (fxx1, fxx, atol=1.0e-10_PREC), &
        'Test LOGADDEXP(y,x) = LOG(EXP(Y) + EXP(X)) for non-critical values')

end subroutine


end
