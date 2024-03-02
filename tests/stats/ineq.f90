
program test_numfort_stats_ineq
    !*  Module with unit tests for inequality measures in stats module.

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_common_testing
    use numfort_stats

    use fcore_testing
    use fcore_strings

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ('stats::ineq unit tests')

    call test_gini (tests)

    call tests%print ()

end subroutine


subroutine test_gini (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(10) :: x, pmf
    real (PREC) :: g, B, g_target
    integer :: i

    tc => tests%add_test ('Gini coeffcient')

    ! Degenerate distribution (this will not work if x = )
    x = 0.0
    pmf = 1
    call gini (x(:1), pmf(:1), g, assume_sorted=.true.)
    call tc%assert_true (g == 0.0, 'Degenerate distribution: Gini = 0')

    ! 2-point distribution where one half has nothing
    x(1:2) = [0.0, 1.0]
    pmf = 0.5
    call gini (x(:2), pmf(:2), g, assume_sorted=.true.)
    call tc%assert_true (abs(g - 0.5_PREC) < 1.0e-10, 'Two-point distribution: Gini = 0.5')

    ! 10-point distribution with same logic as before
    x(1:5) = 0.0
    x(6:10) = 0.1
    pmf = 1.0_PREC / size(pmf)
    call gini (x, pmf, g, assume_sorted=.true.)
    call tc%assert_true (abs(g - 0.5_PREC) < 1.0e-10, '10-point distribution: Gini = 0.5')

    ! Equal distribution among the top N, rest has nothing
    pmf = 1.0_PREC / size(pmf)

    do i = 0, size(x) - 1
        x(1:i) = 0
        x(i+1:) = 1.234
        B = (1.0_PREC - i / real(size(x), PREC)) / 2.0
        g_target = 1.0_PREC - 2.0 * B
        call gini (x, pmf, g, assume_sorted=.true.)
        call tc%assert_true (abs(g-g_target) < 1.0e-10, 'Equal distr. among to Top ' // str (10-i))
    end do


end subroutine


end