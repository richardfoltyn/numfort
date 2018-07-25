

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
    call test_bsearch_cached (tests)

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
    deallocate (haystack)

    n = 1
    allocate (haystack(n), source=0.0_PREC)
    needle = 0.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 1, "size(HAYSTACK) = 1")
    deallocate (haystack)

    n = 11
    allocate (haystack(n))
    call linspace (haystack, 0.0_PREC, 10.0_PREC)

    needle = -1.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1)")

    needle = 11.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK))")

    needle = 0.5d0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 1, "needle in first interval")

    needle = 9.5d0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 10, "needle in last interval")

    deallocate (haystack)

    ! Test with piecewise constant haystack
    allocate (haystack(11))
    call linspace (haystack, 0.0_PREC, 10.0_PREC)
    haystack(3:5) = haystack(3)
    haystack(7:9) = haystack(7)

    needle = 2.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 5, "Piecewise-constant HAYSTACK, lower bound")

    needle = 3.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 5, "Piecewise-constant HAYSTACK, interior value")

    needle = 5.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 6, "Piecewise-constant HAYSTACK, upper bound")

    deallocate (haystack)

end subroutine


subroutine test_bsearch_cached (tests)
    !*  Unit tests for BSEARCH routine
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    real (PREC), dimension(:), allocatable :: haystack
    real (PREC) :: needle
    integer :: i, n
    type (search_cache) :: cache

    tc => tests%add_test ("BSEARCH_CACHED unit tests")

    n = 0
    allocate (haystack(n), source=0.0_PREC)
    needle = 0.0
    call bsearch_cached (needle, haystack, i)
    call tc%assert_true (i == 0, "size(HAYSTACK) = 0")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 0, "size(HAYSTACK) = 0 with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 0, "size(HAYSTACK) = 0 with cache, repeat call")

    deallocate (haystack)

    n = 1
    allocate (haystack(n), source=0.0_PREC)
    needle = 0.0

    call bsearch_cached (needle, haystack, i)
    call tc%assert_true (i == 1, "size(HAYSTACK) = 1")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "size(HAYSTACK) = 1 with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "size(HAYSTACK) = 1 with cache, repeat call")

    deallocate (haystack)

    n = 11
    allocate (haystack(n))
    call linspace (haystack, 0.0_PREC, 10.0_PREC)

    needle = -1.0
    call bsearch_cached (needle, haystack, i)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1)")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1) with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1) with cache, repeat call")

    needle = 11.0
    call bsearch_cached (needle, haystack, i)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK))")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK)) with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK)) with cache, repeat call")

    needle = 0.5d0
    call bsearch_cached (needle, haystack, i)
    call tc%assert_true (i == 1, "needle in first interval")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "needle in first interval with cache")
    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "needle in first interval with cache, repeat call")

    needle = 9.5d0
    call bsearch_cached (needle, haystack, i)
    call tc%assert_true (i == 10, "needle in last interval")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 10, "needle in last interval with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 10, "needle in last interval with cache, repeat call")

    deallocate (haystack)

end subroutine

end program