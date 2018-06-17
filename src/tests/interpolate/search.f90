

program test_numfort_interpolate_search
    !*  Unit tests for search routines of interpolation module

    use, intrinsic :: iso_fortran_env
    use numfort_arrays
    use numfort_common_testing
    use numfort_interpolate

    use fcore_testing

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("Interpolation SEARCH routines unit tests")

    call test_bsearch (tests)

    call tests%print ()

end subroutine


subroutine test_bsearch (tests)
    !*  Unit tests for BSEARCH routine
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(:), allocatable :: haystack
    real (PREC) :: needle
    integer :: i, n

    tc => tests%add_test ("BSEARCH unit tests")

    n = 0
    allocate (haystack(n), source=0.0_PREC)
    needle = 0.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 0, "size(HAYSTACK) = 0")
    i = bsearch (needle, haystack, i)
    call tc%assert_true (i == 0, "size(HAYSTACK) = 0, using init guess")

    deallocate (haystack)

    n = 1
    allocate (haystack(n), source=0.0_PREC)
    needle = 0.0

    i = bsearch (needle, haystack)
    call tc%assert_true (i == 1, "size(HAYSTACK) = 1")
    i = bsearch (needle, haystack, i)
    call tc%assert_true (i == 1, "size(HAYSTACK) = 1, using init guess")

    deallocate (haystack)

    n = 11
    allocate (haystack(n))
    call linspace (haystack, 0.0_PREC, 10.0_PREC)

    needle = -1.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1)")
    i = bsearch (needle, haystack, i)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1), using correct init guess")
    i = bsearch (needle, haystack, i=5)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1), using incorrect init guess")

    needle = 11.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK))")
    i = bsearch (needle, haystack, i)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK)), using correct init guess")
    i = bsearch (needle, haystack, i=4)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK)), using incorrect init guess")

    needle = 0.5d0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 1, "needle in first interval")
    i = bsearch (needle, haystack, i)
    call tc%assert_true (i == 1, "needle in first interval, using correct init guess")
    i = bsearch (needle, haystack, i=3)
    call tc%assert_true (i == 1, "needle in first interval, using incorrect init guess")

    needle = 9.5d0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 10, "needle in last interval")
    i = bsearch (needle, haystack, i)
    call tc%assert_true (i == 10, "needle in last interval, using correct init guess")
    i = bsearch (needle, haystack, i=5)
    call tc%assert_true (i == 10, "needle in last interval, using incorrect init guess")

    ! Test with init guess that is not the correct interval
    needle = 1.5d0
    i = bsearch (needle, haystack, i=5)
    call tc%assert_true (i == 2, "Init guess /= correct bracket")

    deallocate (haystack)

end subroutine

end program