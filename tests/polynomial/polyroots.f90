

program test_polynomial_polyroots
    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_common_testing
    use numfort_polynomial

    use fcore_testing, only: test_suite, test_case
    use fcore_strings

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains


subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("Polynomial root-finding unit tests")

    call test_roots_quad (tests)

    call tests%print ()

end subroutine


subroutine test_roots_quad (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), dimension(3) :: coefs
    real (PREC), dimension(2) :: roots, roots2
    integer :: n
    type (status_t) :: status

    tc => tests%add_test ("POLYROOTS (quad) unit tests")

    ! Try with 2 distinct real roots
    roots = [-1.234d0, 0.234d0]
    coefs(1) = product(roots)
    coefs(2) = - sum(roots)
    coefs(3) = 1.0

    call polyroots (coefs, roots2, n, status)
    call tc%assert_true (all_close (roots, roots2) .and. n == 2, &
        "Polynomial with two distinct real roots")

    ! one unique real root
    roots = 1.2345d0
    coefs(1) = product(roots)
    coefs(2) = - sum(roots)
    coefs(3) = 1.0

    call polyroots (coefs, roots2, n, status)
    call tc%assert_true (all_close (roots, roots2) .and. n == 1, &
        "Polynomial with root with multiplicity two")

    ! No real roots
    coefs(1) = 123.2d0
    coefs(2) = 0.123d0
    coefs(3) = 12.23d0

    call polyroots (coefs, roots2, n, status)
    call tc%assert_true (n==0, "Polynomial with zero real roots")

end subroutine

end program
