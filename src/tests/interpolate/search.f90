

program test_numfort_interpolate_search
    !*  Unit tests for search routines of interpolation module

    use, intrinsic :: iso_fortran_env
    use numfort_arrays
    use numfort_common_testing
    use numfort_interpolate
    use numfort_stats

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

    call test_interp_find (tests)

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
    call tc%assert_true (i == -1, "size(HAYSTACK) = 0")
    deallocate (haystack)

    n = 1
    allocate (haystack(n), source=0.0_PREC)
    needle = 0.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == -1, "size(HAYSTACK) = 1")
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
    i = bsearch (needle, haystack)
    call tc%assert_true (i == -1, "size(HAYSTACK) = 0")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == -1, "size(HAYSTACK) = 0 with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == -1, "size(HAYSTACK) = 0 with cache, repeat call")

    deallocate (haystack)

    n = 1
    allocate (haystack(n), source=0.0_PREC)
    needle = 0.0

    i = bsearch (needle, haystack)
    call tc%assert_true (i == -1, "size(HAYSTACK) = 1")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == -1, "size(HAYSTACK) = 1 with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == -1, "size(HAYSTACK) = 1 with cache, repeat call")

    deallocate (haystack)

    n = 11
    allocate (haystack(n))
    call linspace (haystack, 0.0_PREC, 10.0_PREC)

    needle = -1.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1)")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1) with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "needle < HAYSTACK(1) with cache, repeat call")

    needle = 11.0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK))")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK)) with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 10, "needle > HAYSTACK(size(HAYSTACK)) with cache, repeat call")

    needle = 0.5d0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 1, "needle in first interval")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "needle in first interval with cache")
    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 1, "needle in first interval with cache, repeat call")

    needle = 9.5d0
    i = bsearch (needle, haystack)
    call tc%assert_true (i == 10, "needle in last interval")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 10, "needle in last interval with cache")

    call bsearch_cached (needle, haystack, i, cache)
    call tc%assert_true (i == 10, "needle in last interval with cache, repeat call")

    deallocate (haystack)

end subroutine



subroutine test_interp_find (tests)
    !*  Unit tests for INTERP_FIND routine
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (search_cache) :: cache
    type (status_t) :: status
    real (PREC), dimension(:), allocatable :: haystack, needles, rwork
    real (PREC), dimension(:), allocatable :: weight, weight_ok
    integer, dimension(:), allocatable :: ilbound, iorder, ilbound_ok
    integer :: i, ilb, n
    real (PREC) :: xi, wgt
    logical :: all_ok

    tc => tests%add_test ("INTERP_FIND unit tests")

    ! === Input validation ===

    ! --- Invalid knots array size ---

    allocate (needles(2), ilbound(2), weight(2))

    do i = 0, 1
        allocate (haystack(i), source=0.0_PREC)
        ! Test scalar interface
        status = NF_STATUS_UNDEFINED
        call interp_find (xi, haystack, ilb, wgt, status=status)
        call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
            "0d API: Invalid KNOTS argument")

        ! Test 1d interface
        status = NF_STATUS_UNDEFINED
        call interp_find (needles, haystack, ilbound, weight, status=status)
        call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
            "1d API: Invalid KNOTS argument")

        deallocate (haystack)
    end do

    deallocate (needles, ilbound, weight)

    ! --- Incompatible array sizes ---

    allocate (needles(1), ilbound(2), weight(2))
    allocate (haystack(5), source=0.0_PREC)

    status = NF_STATUS_UNDEFINED
    call interp_find (needles, haystack, ilbound, weight, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable NEEDLES, ILBOUND arguments")

    deallocate (needles, ilbound, weight, haystack)

    ! Non-conformable ILBOUND, WEIGHT arrays
    allocate (needles(2), ilbound(2), weight(1))
    allocate (haystack(10))

    status = NF_STATUS_UNDEFINED
    call interp_find (needles, haystack, ilbound, weight, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Non-conformable ILBOUND, WEIGHT arguments")

    deallocate (needles, ilbound, weight, haystack)

    ! === identical Needles and Knots ===
    n = 11
    allocate (haystack(n), rwork(n), ilbound(n), ilbound_ok(n), weight(n))
    allocate (iorder(n))

    call set_seed (1234)

    call random_number (rwork)
    call argsort (rwork, iorder)
    forall (i=1:n) haystack(i) = rwork(iorder(i))

    allocate (needles(n), source=haystack)

    ! Correct indices
    call arange (ilbound_ok, 1)
    ilbound_ok(n) = n-1

    ! Scalar interface
    do i = 1, size(needles)
        xi = needles(i)
        status = NF_STATUS_UNDEFINED
        call interp_find (xi, haystack, ilbound(i), weight(i), status=status)
    end do

    all_ok = all(ilbound == ilbound_ok) .and. all(weight(1:n-1) == 1.0_PREC)
    all_ok = all_ok .and. (weight(n) == 0.0_PREC)
    call tc%assert_true (all_ok .and. status == NF_STATUS_OK, &
        '0d API: NEEDLES identical to KNOTS')

    ! Array interface
    ilbound(:) = 0
    weight(:) = 0.0
    status = NF_STATUS_UNDEFINED
    call interp_find (needles, haystack, ilbound, weight, status=status)

    all_ok = all(ilbound == ilbound_ok) .and. all(weight(1:n-1) == 1.0_PREC)
    all_ok = all_ok .and. (weight(n) == 0.0_PREC)

    deallocate (rwork, iorder)
    deallocate (needles, haystack, ilbound, ilbound_ok, weight)

    ! === Interpolating interior/non-interior points ===
    n = 10
    allocate (needles(n+1), haystack(n), ilbound(n+1), ilbound_ok(n+1))
    allocate (weight(n+1), weight_ok(n+1))

    call arange(haystack, 1.0_PREC)
    needles(1) = 0.5_PREC
    needles(2:n) = (haystack(2:n) + haystack(1:n-1)) / 2.0_PREC
    needles(n+1) = 10.5
    ! Throw in some arbitrary non-interior points
    needles(4) = 0.0
    needles(7) = 11.0

    call arange(ilbound_ok(2:n-1), 1)
    ilbound_ok(1) = 1
    ilbound_ok(n:n+1) = n-1
    ilbound_ok(4) = 1
    ilbound_ok(7) = n-1

    weight_ok(:) = 0.5_PREC
    weight_ok(1) = 1.5_PREC
    weight_ok(n+1) = -0.5_PREC
    weight_ok(4) = 2.0_PREC
    weight_ok(7) = -1.0_PREC

    ! Scalar interface, default extrapolation
    do i = 1, size(needles)
        status = NF_STATUS_UNDEFINED
        xi = needles(i)
        call interp_find (xi, haystack, ilbound(i), weight(i), status=status)
    end do

    all_ok = all(ilbound == ilbound_ok)
    all_ok = all_ok .and. all_close (weight, weight_ok, rtol=0.0_PREC, atol=1.0e-12_PREC)

    call tc%assert_true (all_ok .and. status == NF_STATUS_OK, &
        "0d API: Mix of interior and non-interior interpolation points")

    ! 1d interface, default extrapolation
    status = NF_STATUS_UNDEFINED
    weight(:) = 0.0
    ilbound(:) = 0

    call interp_find (needles, haystack, ilbound, weight, status=status)

    all_ok = all(ilbound == ilbound_ok)
    all_ok = all_ok .and. all_close (weight, weight_ok, rtol=0.0_PREC, atol=1.0e-12_PREC)

    call tc%assert_true (all_ok .and. status == NF_STATUS_OK, &
        "1d API: Mix of interior and non-interior interpolation points ")

    ! === No extrapolation, relace with boundary values ===

    where (ilbound_ok < 1)
        ilbound_ok = 1
    end where

    where (weight_ok > 1.0)
        weight_ok = 1.0
    else where (weight_ok < 0.0)
        weight_ok = 0.0
    end where

    ! Scalar interface
    do i = 1, size(needles)
        status = NF_STATUS_UNDEFINED
        xi = needles(i)
        call interp_find (xi, haystack, ilbound(i), weight(i), &
            ext=NF_INTERP_EVAL_BOUNDARY, status=status)
    end do

    all_ok = all(ilbound == ilbound_ok)
    all_ok = all_ok .and. all_close (weight, weight_ok, rtol=0.0_PREC, atol=1.0e-12_PREC)

    call tc%assert_true (all_ok .and. status == NF_STATUS_OK, &
        "0d API: Some non-interior points, ext=NF_INTERP_EVAL_BOUNDARY")

    ! 1d interface
    status = NF_STATUS_UNDEFINED
    weight(:) = 0.0
    ilbound(:) = 0

    call interp_find (needles, haystack, ilbound, weight, &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)

    all_ok = all(ilbound == ilbound_ok)
    all_ok = all_ok .and. all_close (weight, weight_ok, rtol=0.0_PREC, atol=1.0e-12_PREC)

    call tc%assert_true (all_ok .and. status == NF_STATUS_OK, &
        "1d API: Some non-interior points, ext=NF_INTERP_EVAL_BOUNDARY ")

    ! === Error on non-interior point ===

    all_ok = .true.

    ! Scalar interface
    do i = 1, size(needles)
        status = NF_STATUS_UNDEFINED
        xi = needles(i)
        call interp_find (xi, haystack, ilbound(i), weight(i), &
            ext=NF_INTERP_EVAL_ERROR, status=status)
        if (xi < haystack(1) .or. xi > haystack(n)) then
            all_ok = all_ok .and. (status == NF_STATUS_BOUNDS_ERROR)
        else
            all_ok = all_ok .and. (status == NF_STATUS_OK)
        end if
    end do

    call tc%assert_true (all_ok, "0d API: Some non-interior points, ext=NF_INTERP_EVAL_ERROR")

    ! 1d interface
    status = NF_STATUS_UNDEFINED

    call interp_find (needles, haystack, ilbound, weight, &
        ext=NF_INTERP_EVAL_ERROR, status=status)

    call tc%assert_true (status == NF_STATUS_BOUNDS_ERROR, &
        "1d API: Some non-interior points, ext=NF_INTERP_EVAL_ERROR")

    ! === Replace non-interior points with const ===
    ! We assume that in this case the FIND routine returns specific (invalid)
    ! lower bound indices to signal to the EVAL routine that these should
    ! be replaced with some left/right constants

    where (needles < haystack(1))
        ilbound_ok = 0
    else where (needles > haystack(n))
        ilbound_ok = n
    end where

    ! Scalar interface
    ilbound(:) = -1
    weight(:) = 0.0

    all_ok = .true.
    do i = 1, size(needles)
        status = NF_STATUS_UNDEFINED
        xi = needles(i)
        call interp_find (xi, haystack, ilbound(i), weight(i), &
            ext=NF_INTERP_EVAL_CONST, status=status)
        all_ok = all_ok .and. (status == NF_STATUS_OK)
    end do

    all_ok = all_ok .and. all(ilbound == ilbound_ok)
    call tc%assert_true (all_ok, &
        "0d API: Some non-interior points, ext=NF_INTERP_EVAL_CONST")

    status = NF_STATUS_UNDEFINED
    call interp_find (needles, haystack, ilbound, weight, &
        ext=NF_INTERP_EVAL_CONST, status=status)

    all_ok = all(ilbound == ilbound_ok) .and. (status == NF_STATUS_OK)

    call tc%assert_true (all_ok, &
        "1d API: Some non-interior points, ext=NF_INTERP_EVAL_CONST")

end subroutine

end program
