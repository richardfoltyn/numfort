


pure function __APPEND(bsearch,__PREC) (needle, haystack) result(ilb)
    !*  BSEARCH finds the interval on array HAYSTACK that contains NEEDLE,
    !   returning the array index of the lower boundary of that interval.
    !
    !   If NEEDLE is smaller than the first element of HAYSTACK, the rountine
    !   returns 1. If NEEDLE is larger than the last of element of HAYSTACK,
    !   size(HAYSTACK)-1 is returned.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: needle
        !*  Value to bracket
    real (PREC), intent(in), dimension(:), contiguous :: haystack
        !*  Array of sorted (!) values that are interpreted as the boundaries
        !   of a sequence of bracketing intervals.
    integer :: ilb
        !*  Index of array element that is the lower bound of the bracketing
        !   interval such that HAYSTACK(ilb) <= NEEDLE < HAYSTACK(ilb+1),

    integer :: iub, imid, n

    n = size(haystack)
    if (n <= 1) then
        ilb = -1
        return
    end if

    ilb = 1
    iub = n

    do while (iub > (ilb + 1))
        imid = (iub + ilb) / 2
        if (haystack(imid) > needle) then
            iub = imid
        else
            ilb = imid
        end if
    end do

end function



pure subroutine __APPEND(bsearch_cached,__PREC) (needle, haystack, ilb, cache)
    !*  BSEARCH finds the interval on array HAYSTACK that contains NEEDLE,
    !   returning the array index of the lower boundary of that interval.
    !
    !   If NEEDLE is smaller than the first element of HAYSTACK, the rountine
    !   returns 1. If NEEDLE is larger than the last of element of HAYSTACK,
    !   size(HAYSTACK)-1 is returned.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: needle
        !*  Value to bracket
    real (PREC), intent(in), dimension(:), contiguous :: haystack
        !*  Array of sorted (!) values that are interpreted as the boundaries
        !   of a sequence of bracketing intervals.
    integer, intent(out) :: ilb
        !*  Index of array element that is the lower bound of the bracketing
        !   interval such that HAYSTACK(ilb) <= NEEDLE < HAYSTACK(ilb+1),
        !   if such a bracket exists.
    type (search_cache), intent(inout) :: cache
        !*  Optional search cache. If present, initially the routine checks
        !   whether the interval HAYSTACK(i) <= NEEDLE < HAYSTACK(i+1)),
        !   where i is the cached index and returns immediately
        !   if this condition holds. Otherwise the regular search algorithm
        !   is resumed.

    integer :: iub, imid, n, i

    n = size(haystack)
    if (n <= 1) then
        ilb = -1
        return
    end if

    i = cache%i
    i = max(min(i, n-1), 1)
    if (haystack(i) <= needle) then
        ilb = i
        if (haystack(i+1) > needle) then
            goto 100
        else if (i == (n-1)) then
            goto 100
        end if
        iub = n
    else
        ilb = 1
        iub = i
    end if

    do while (iub > (ilb + 1))
        imid = (iub + ilb) / 2
        if (haystack(imid) > needle) then
            iub = imid
        else
            ilb = imid
        end if
    end do

100 continue
    cache%i = ilb

end subroutine


pure function __APPEND(interp_find,__PREC) (needle, haystack) result (res)
    !*  INTERP_FIND is a wrapper around BSEARCH that exists only for
    !   backwards compatibility.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:), contiguous :: haystack
    real (PREC), intent(in) :: needle

    integer :: res

    res = bsearch (needle, haystack)

end function



pure subroutine __APPEND(interp_find_cached,__PREC) (needle, haystack, ilb, cache)
    !*  INTERP_FIND_CACHED is a wrapper around BSEARCH that exists only for
    !   backwards compatibility.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:), contiguous :: haystack
    integer, intent(out) :: ilb
    type (search_cache), intent(inout) :: cache

    call bsearch_cached (needle, haystack, ilb, cache)

end subroutine


pure subroutine __APPEND(interp_find_wgt_cached,__PREC) (needle, haystack, &
        ilb, wgt_lb, cache)
    !*  INTERP_FIND_WGT_CACHED returns the index and weight associated with
    !   the lower bound of a (bracketing) interpolation interval, optionally
    !   using a search cache.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:), contiguous :: haystack
    integer, intent(out) :: ilb
        !*  Index of the lower bound of the interpolation bracket
    real (PREC), intent(out) :: wgt_lb
        !*  Interpolation weight on the lower bound of the interpolation
        !   bracket.
    type (search_cache), intent(inout) :: cache

    call interp_find_cached (needle, haystack, ilb, cache)

    ! Compute interpolation weight
    wgt_lb = (haystack(ilb+1)-needle) / (haystack(ilb+1)-haystack(ilb))
end subroutine


pure subroutine __APPEND(interp_find_wgt_cached_1d,__PREC) (needle, haystack, &
       ilb, wgt_lb)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:), contiguous :: needle
    real (PREC), intent(in), dimension(:), contiguous :: haystack
    integer, intent(out), dimension(:), contiguous :: ilb
        !*  Index of the lower bound of the interpolation bracket
    real (PREC), intent(out), dimension(:), contiguous :: wgt_lb
        !*  Interpolation weight on the lower bound of the interpolation
        !   bracket.

    type (search_cache) :: cache
    real (PREC) :: x
    integer :: i, j

    do i = 1, size(needle)
        x = needle(i)
        call bsearch_cached (x, haystack, j, cache)
        wgt_lb(i) = (haystack(j+1)-x) / (haystack(j+1)-haystack(j))
        ilb(i) = j
    end do
end subroutine
