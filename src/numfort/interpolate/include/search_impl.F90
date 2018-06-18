


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
    real (PREC), intent(in), dimension(:) :: haystack
        !*  Array of sorted (!) values that are interpreted as the boundaries
        !   of a sequence of bracketing intervals.
    integer :: ilb
        !*  Index of array element that is the lower bound of the bracketing
        !   interval such that HAYSTACK(ilb) <= NEEDLE < HAYSTACK(ilb+1),

    integer :: iub, imid, n

    n = size(haystack)

    ilb = min(1, n)
    if (n <= 1) return

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
    real (PREC), intent(in), dimension(:) :: haystack
        !*  Array of sorted (!) values that are interpreted as the boundaries
        !   of a sequence of bracketing intervals.
    integer, intent(out) :: ilb
        !*  Index of array element that is the lower bound of the bracketing
        !   interval such that HAYSTACK(ilb) <= NEEDLE < HAYSTACK(ilb+1),
        !   if such a bracket exists.
    type (search_cache), intent(inout), optional :: cache
        !*  Optional search cache. If present, initially the routine checks
        !   whether the interval HAYSTACK(i) <= NEEDLE < HAYSTACK(i+1)),
        !   where i is the cached index and returns immediately
        !   if this condition holds. Otherwise the regular search algorithm
        !   is resumed.

    integer :: iub, imid, n, i

    n = size(haystack)

    ilb = min(1, n)
    if (n <= 1) goto 100

    if (present(cache)) then
        i = cache%i
        i = max(min(i, n), 1)
        if (haystack(i) <= needle) then
            ilb = i
            if (i == (n-1)) then
                goto 100
            else if (haystack(i+1) > needle) then
                goto 100
            end if
            iub = n
        else
            ilb = 1
            iub = i
        end if
    else
        ilb = 1
        iub = n
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
    if (present(cache)) cache%i = ilb

end subroutine


pure function __APPEND(interp_find,__PREC) (needle, haystack) result (res)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: haystack
    real (PREC), intent(in) :: needle

    integer :: res, n

    n = size(haystack)

    if (n < 2) then
        res = -1
    else if (needle <= haystack(1)) then
        res = 1
    else if (needle >= haystack(n-1)) then
        res = n-1
    else
        res = bsearch (needle, haystack)
    end if

end function



pure subroutine __APPEND(interp_find_cached,__PREC) (needle, haystack, ilb, cache)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: haystack
    real (PREC), intent(in) :: needle
    integer, intent(out) :: ilb
    type (search_cache), intent(inout), optional :: cache

    integer :: n

    n = size(haystack)

    if (n < 2) then
        ilb = -1
    else if (needle <= haystack(1)) then
        ilb = 1
    else if (needle >= haystack(n-1)) then
        ilb = n-1
    else
        call bsearch_cached (needle, haystack, ilb, cache)
        ! Cache is updated in BSEARCH_CACHED, we can skip it here
        return
    end if

    if (present(cache)) cache%i = ilb

end subroutine

