

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
    call test_interp_find_decr (tests)

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




subroutine test_interp_find_decr (tests)
    !*  Unit tests for INTERP_FIND_DECR which performs searches on
    !   (weakly) decreasing input data.
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    type (status_t) :: status
    real (PREC), dimension(:), allocatable :: x, weight, knots, rwork, weight_ok
    integer, dimension(:), allocatable :: ilbound, ilbound_ok
    integer :: nx, nk
    logical :: values_ok

    tc => tests%add_test ('Unit tests for INTERP_FIND_DECR')

    call set_seed (12345)

    ! === Input checks ===

    ! 1. Invalid number of knots
    nx = 10
    nk = 1

    allocate (x(nx), weight(nx), ilbound(nx), knots(nk))
    call linspace (x, 1.0_PREC, 10.0_PREC)

    status = NF_STATUS_UNDEFINED

    call interp_find_decr (x, knots, ilbound, weight, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'Invalid HAYSTACK argument')

    deallocate (x, weight, ilbound, knots)

    ! 2. Non-conformable X, ILBOUND
    nx = 3
    nk = 2
    allocate (x(nx), ilbound(nx+1), weight(nx+1), knots(nk))
    call linspace (x, 10.0_PREC, 1.0_PREC)

    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'Non-conformable NEEDLE, ILBOUND')

    deallocate (x, ilbound, weight, knots)

    ! 2. Non-conformable ILBOUND, WEIGHT
    nx = 5
    nk = 11
    allocate (x(nx), ilbound(nx), weight(nx-1), knots(nk))
    call linspace (x, 10.0_PREC, 1.0_PREC)

    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'Non-conformable ILBOUND, WEIGHT')

    deallocate (x, ilbound, weight, knots)

    ! 3. Increasing values in X
    nx = 11
    nk = 5
    allocate (x(nx), ilbound(nx), weight(nx), knots(nk))

    ! Test with interior violation
    call linspace (x, 10.0_PREC, 1.0_PREC)
    x(5) = 50.0
    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'Increasing value in X, #1')

    ! Test with violation at the beginning
    call linspace (x, 10.0_PREC, 1.0_PREC)
    x(1) = 1.0
    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'Increasing value in X, #2')

    ! Test with violation at end
    call linspace (x, 10.0_PREC, 1.0_PREC)
    x(nx) = 10.0
    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, 'Increasing value in X, #3')

    deallocate (x, ilbound, weight, knots)

    ! === Evaluate with valid input data ===
    nx = 17
    nk = 11
    allocate (x(nx), ilbound(nx), weight(nx), knots(nk))
    allocate (ilbound_ok(nx), weight_ok(nx))
    call linspace (knots, 1.0_PREC, 10.0_PREC)
    allocate (rwork(nk))
    call random_number (rwork)
    ! Add some perturbation that does not alter ordering in KNOTS
    knots(:) = knots + rwork * 0.5

    ! === Evaluate with all interior points ===
    call linspace (x, 10.0_PREC, 1.0_PREC)

    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status)
    ! Compute with standard routine
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and. &
        all_close (weight, weight_ok, rtol=0.0_PREC, atol=1.0e-14_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Evaluating with all-interior points')

    ! === Evaluate with points outside of KNOTS ===
    call linspace (x, 11.0_PREC, -1.0_PREC)

    status = NF_STATUS_UNDEFINED
    ilbound(:) = 0
    weight(:) = 0.0
    call interp_find_decr (x, knots, ilbound, weight, status)
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and. &
        all_close (weight, weight_ok, atol=1.0e-14_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Evaluating with interior/exterior points')

    deallocate (x, ilbound, ilbound_ok, weight, weight_ok, knots, rwork)

    ! === Long KNOTS argument ===
    ! Test with longs knots argument such that this triggers binary search
    nk = 101
    nx = 10

    allocate (knots(nk), rwork(nk))
    allocate (x(nx), ilbound(nx), ilbound_ok(nx), weight(nx), weight_ok(nx))

    call linspace (knots, 0.0_PREC, 100.0_PREC)
    call random_number (rwork)
    knots(:) = knots + 0.5 * rwork
    deallocate (rwork)

    ! 1. All interior points
    call linspace (x, 99.0_PREC, 1.0_PREC)
    allocate (rwork(nx))
    call random_number (rwork)
    x(:) = x + rwork

    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status)
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and.  &
        all_close (weight, weight_ok, atol=1.0e-14_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Long HAYSTACK input, all interior X')

    ! 2. Exterior/interior X
    call linspace (x, 200.0_PREC, -10.0_PREC)
    call random_number (rwork)
    x(:) = x + 5.0 * rwork

    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status)
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and.  &
        all_close (weight, weight_ok, atol=1.0e-14_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Long HAYSTACK input, all exterior X')

    deallocate (x, knots, rwork, ilbound, ilbound_ok, weight, weight_ok)

    ! === Constant function ===

    nk = 11
    nx = 17
    allocate (x(nx), ilbound(nx), ilbound_ok(nx), weight(nx), weight_ok(nx))
    allocate (knots(nk))

    call linspace (knots, -10.0_PREC, 10.0_PREC)
    allocate (rwork(nk))
    call random_number (rwork)
    knots(:) = knots + rwork * 0.1
    deallocate (rwork)

    ! Interior within HAYSTACK
    x(1:3) = knots(5)
    x(4:) = knots(4) - 0.1
    x(7:) = knots(3) - 0.1

    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status)
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and. &
        all_close (weight, weight_ok, atol=1.0e-14_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Constant X at interior of HAYSTACK')

    ! Left bound of HAYSTACK
    x(1) = knots(1)
    x(2:4) = knots(1) - 1.0
    x(5:) = knots(1) - 2.0
    status = NF_STATUS_UNDEFINED
    ilbound(:) = 0
    weight(:) = -1.0
    call interp_find_decr (x, knots, ilbound, weight, status)
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and. &
        all_close (weight, weight_ok, atol=1.0e-14_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Constant X, left of HAYSTACK(1)')

    ! Right of HAYSTACK
    x(1) = knots(nk) + 1.0
    x(2:5) = knots(nk) + 0.5
    x(6:) = knots(nk) - 1.0
    status = NF_STATUS_UNDEFINED
    ilbound(:) = 0
    weight(:) = -1.0
    call interp_find_decr (x, knots, ilbound, weight, status)
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and. &
        all_close (weight, weight_ok, atol=1.0e-14_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Constant X, right of HAYSTACK(-1)')

    deallocate (x, ilbound, ilbound_ok, weight, weight_ok, knots)

    ! === Constant function, large HAYSTACK ===

    nk = 101
    nx = 11
    allocate (x(nx), ilbound(nx), ilbound_ok(nx), weight(nx), weight_ok(nx))
    allocate (knots(nk))

    call linspace (knots, -10.0_PREC, 10.0_PREC)
    allocate (rwork(nk))
    call random_number (rwork)
    knots(:) = knots + rwork * 0.05
    deallocate (rwork)

    ! Interior within HAYSTACK
    x(1:3) = knots(5)
    x(4:) = knots(4) - 0.1
    x(7:) = knots(3) - 0.1

    status = NF_STATUS_UNDEFINED
    call interp_find_decr (x, knots, ilbound, weight, status)
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and. &
        all_close (weight, weight_ok, atol=1.0e-14_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Constant X at interior of HAYSTACK')

    ! Left bound of HAYSTACK
    x(1) = knots(1)
    x(2:4) = knots(1) - 1.0
    x(5:) = knots(1) - 2.0
    status = NF_STATUS_UNDEFINED
    ilbound(:) = 0
    weight(:) = -1.0
    call interp_find_decr (x, knots, ilbound, weight, status)
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and. &
        all_close (weight, weight_ok, atol=1.0e-14_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Constant X, left of HAYSTACK(1)')

    ! Right of HAYSTACK
    x(1) = knots(nk) + 1.0
    x(2:5) = knots(nk) + 0.5
    x(6:) = knots(nk) - 1.0
    status = NF_STATUS_UNDEFINED
    ilbound(:) = 0
    weight(:) = -1.0
    call interp_find_decr (x, knots, ilbound, weight, status)
    call interp_find (x, knots, ilbound_ok, weight_ok)
    values_ok = all(ilbound == ilbound_ok) .and. &
        all_close (weight, weight_ok, atol=1.0e-14_PREC, rtol=0.0_PREC)
    call tc%assert_true (status == NF_STATUS_OK .and. values_ok, &
        'Constant X, right of HAYSTACK(-1)')

    deallocate (x, ilbound, ilbound_ok, weight, weight_ok, knots)



end subroutine

end program
