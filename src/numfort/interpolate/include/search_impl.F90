


pure subroutine bsearch (needle, haystack, ilb)
    !*  BSEARCH finds the interval on array HAYSTACK that contains NEEDLE,
    !   returning the array index of the lower boundary of that interval.
    !
    !   If NEEDLE is smaller than the first element of HAYSTACK, the rountine
    !   returns 1. If NEEDLE is larger than the last of element of HAYSTACK,
    !   size(HAYSTACK)-1 is returned.
    real (PREC), intent(in) :: needle
        !*  Value to bracket
    real (PREC), intent(in), dimension(:), contiguous :: haystack
        !*  Array of sorted (!) values that are interpreted as the boundaries
        !   of a sequence of bracketing intervals.
    integer, intent(out) :: ilb
        !*  Index of array element that is the lower bound of the bracketing
        !   interval such that HAYSTACK(ilb) <= NEEDLE < HAYSTACK(ilb+1),

    if (size(haystack) <= 1) then
        ilb = -1
        return
    end if

    call bsearch_impl (needle, haystack, ilb)

end subroutine



pure subroutine bsearch_impl (needle, haystack, ilb)
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:), contiguous :: haystack
    integer, intent(out) :: ilb

    integer :: iub, imid


    ilb = 1
    iub = size(haystack)

    do while (iub > (ilb + 1))
        imid = (iub + ilb) / 2
        if (haystack(imid) > needle) then
            iub = imid
        else
            ilb = imid
        end if
    end do

end subroutine



pure subroutine bsearch_1d (needle, haystack, ilbound)
    !*  BSEARCH finds, for each element of NEEDLE, the interval on array
    !   HAYSTACK that contains the element, returning the array index of
    !   the lower boundary of that interval.
    !
    !   If a value in NEEDLE is smaller than the first element of HAYSTACK,
    !   the rountine returns 1 for that element.
    !   If NEEDLE is larger than the last of element of HAYSTACK,
    !   size(HAYSTACK)-1 is returned.
    real (PREC), intent(in), dimension(:), contiguous :: needle
        !*  Array of values to bracket
    real (PREC), intent(in), dimension(:), contiguous :: haystack
        !*  Array of sorted (!) values that are interpreted as the boundaries
        !   of a sequence of bracketing intervals.
    integer, intent(out), dimension(:), contiguous :: ilbound

    integer :: i, n

    if (size(haystack) <= 1) then
        ilbound = -1
        return
    end if

    n = size(needle)

    do i = 1, n
        call bsearch_impl (needle(i), haystack, ilbound(i))
    end do


end subroutine



pure subroutine bsearch_cached_1d (needle, haystack, ilbound, cache)
    !*  BSEARCH finds, for each element of NEEDLE, the interval on array
    !   HAYSTACK that contains the element, returning the array index of
    !   the lower boundary of that interval.
    !
    !   If a value in NEEDLE is smaller than the first element of HAYSTACK,
    !   the rountine returns 1 for that element.
    !   If NEEDLE is larger than the last of element of HAYSTACK,
    !   size(HAYSTACK)-1 is returned.
    real (PREC), intent(in), dimension(:), contiguous :: needle
        !*  Array of values to bracket
    real (PREC), intent(in), dimension(:), contiguous :: haystack
        !*  Array of sorted (!) values that are interpreted as the boundaries
        !   of a sequence of bracketing intervals.
    integer, intent(out), dimension(:), contiguous :: ilbound
    type (search_cache), intent(inout) :: cache

    integer :: i, n, ilb

    if (size(haystack) <= 1) then
        ilbound = -1
        return
    end if

    n = size(needle)
    ilb = max(min(cache%i, n-1), 1)

    do i = 1, n
        call bsearch_cached_impl (needle(i), haystack, ilb)
        ilbound(i) = ilb
    end do

    cache%i = ilb

end subroutine



pure subroutine bsearch_cached (needle, haystack, ilb, cache)
    !*  BSEARCH_CACHED finds the interval on array HAYSTACK that contains NEEDLE,
    !   returning the array index of the lower boundary of that interval.
    !
    !   If NEEDLE is smaller than the first element of HAYSTACK, the rountine
    !   returns 1. If NEEDLE is larger than the last of element of HAYSTACK,
    !   size(HAYSTACK)-1 is returned.
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

    integer :: n

    n = size(haystack)
    if (n <= 1) then
        ilb = -1
        return
    end if

    ilb = max(min(cache%i, n-1), 1)

    call bsearch_cached_impl (needle, haystack, ilb)

    cache%i = ilb
end subroutine



pure subroutine bsearch_cached_impl (needle, haystack, ilb)
    !*  BSEARCH_CACHED_IMPL finds the interval on array HAYSTACK that contains
    !   NEEDLE, returning the array index of the lower boundary of that interval.
    !
    !   If NEEDLE is smaller than the first element of HAYSTACK, the rountine
    !   returns 1. If NEEDLE is larger than the last of element of HAYSTACK,
    !   size(HAYSTACK)-1 is returned.
    !
    !   Note: Implementation routine, does not perform any input checks!
    real (PREC), intent(in) :: needle
        !*  Value to bracket
    real (PREC), intent(in), dimension(:), contiguous :: haystack
        !*  Array of sorted (!) values that are interpreted as the boundaries
        !   of a sequence of bracketing intervals.
    integer, intent(inout) :: ilb
        !*  Index of array element that is the lower bound of the bracketing
        !   interval such that HAYSTACK(ilb) <= NEEDLE < HAYSTACK(ilb+1),
        !   if such a bracket exists.

    integer :: iub, imid, n

    if (haystack(ilb) <= needle) then
        n = size(haystack)
        if (haystack(ilb+1) > needle) then
            return
        else if (ilb == (n-1)) then
            return
        end if
        iub = n
    else
        iub = ilb
        ilb = 1
    end if

    do while (iub > (ilb + 1))
        imid = (iub + ilb) / 2
        if (haystack(imid) > needle) then
            iub = imid
        else
            ilb = imid
        end if
    end do

end subroutine



pure function bsearch_bounded (needle, haystack, ilb, iub) result(res)
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:), contiguous :: haystack
    integer, intent(in) :: ilb, iub
    integer :: res

    integer :: imid, lilb, liub

    lilb = ilb
    liub = iub

    do while (liub > (lilb + 1))
        imid = (liub + lilb) / 2
        if (haystack(imid) > needle) then
            liub = imid
        else
            lilb = imid
        end if
    end do

    res = lilb

end function



pure function lsearch_bounded (needle, haystack, ilb, iub) result(res)
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:), contiguous :: haystack
    integer, intent(in) :: ilb, iub
    integer :: res

    res = ilb

    do res = ilb, iub-2
        if (needle < haystack(res+1)) exit
    end do

end function



pure subroutine interp_find_check_input (x, knots, ilbound, weight, status)
    !*  INTERP_FIND_CHECK_INPUT performs input validation for various
    !   INTERP_FIND routines.

    real (PREC), intent(in), dimension(:), optional :: x
        !*  Array of points for which to determine the location within `knots`
    real (PREC), intent(in), dimension(:) :: knots
        !*  Array defining the points over which to search.
    integer, intent(in), dimension(:), optional :: ilbound
        !*  Array to store indices of bracket lower bounds. Optional since
        !   this is not applicable for the scalar API.
    real (PREC), intent(in), dimension(:), optional :: weight
        !*  Array to store weights on bracket lower bounds. Optional since
        !   this is not applicable for the scalar API.
    type (status_t), intent(out) :: status

    status = NF_STATUS_OK

    if (size(knots) < 2) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(x) .and. present(ilbound) .and. present(weight)) then
        if (size(x) /= size(ilbound)) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if

        if (size(ilbound) /= size(weight)) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

100 continue

end subroutine



pure subroutine interp_find_cached_impl_0d (x, knots, ilbound, weight, &
        ext, cache, status)
    !*  INTERP_FIND_IMPL_CACHE_SCALAR implements an algorithm to identify the
    !   bracketing interval and interpolation weights that are required for
    !   various interpolation methods.
    real (PREC), intent(in) :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out) :: ilbound
    real (PREC), intent(out) :: weight
    integer (NF_ENUM_KIND), intent(in) :: ext
    type (search_cache), intent(inout) :: cache
        !*  Cache to speed up search algorithm.
    type (status_t), intent(out) :: status

    integer :: nknots, j

    status = NF_STATUS_OK

    nknots = size(knots)

    if (x < knots(1)) then
        select case (ext)
        case (NF_INTERP_EVAL_BOUNDARY)
            ilbound = 1
            weight = 1.0_PREC
            return
        case (NF_INTERP_EVAL_CONST)
            ilbound = 0
            ! Weight does not matter in this case, should not be used
            weight = 0.0_PREC
            return
        case (NF_INTERP_EVAL_ERROR)
            status = NF_STATUS_BOUNDS_ERROR
            return
        case (NF_INTERP_EVAL_EXTRAPOLATE)
            ilbound = 1
            weight = (knots(2) - x) / (knots(2) - knots(1))
            return
        end select
    else if (x > knots(nknots)) then
        select case (ext)
        case (NF_INTERP_EVAL_BOUNDARY)
            ! Upper bound: set corresponding interpolation interval to the
            ! last one. Correct boundary value will be interpolated below.
            ilbound = nknots - 1
            weight = 0.0_PREC
            return
        case (NF_INTERP_EVAL_CONST)
            ilbound = nknots
            ! Weight does not matter in this case, should not be used
            weight = 0.0_PREC
            return
        case (NF_INTERP_EVAL_ERROR)
            status = NF_STATUS_BOUNDS_ERROR
            return
        case (NF_INTERP_EVAL_EXTRAPOLATE)
            ilbound = nknots - 1
            weight = (knots(nknots) - x) / (knots(nknots) - knots(nknots-1))
            return
        end select
    end if

    call bsearch_cached (x, knots, j, cache)

    ilbound = j
    weight = (knots(j+1) - x) / (knots(j+1) - knots(j))

end subroutine



pure subroutine interp_find_impl_0d (x, knots, ilbound, weight, ext, status)
    !*  INTERP_FIND_IMPL_SCALAR implements an algorithm to identify the
    !   bracketing interval and interpolation weights that are required for
    !   various interpolation methods.
    real (PREC), intent(in) :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out) :: ilbound
    real (PREC), intent(out) :: weight
    integer (NF_ENUM_KIND), intent(in) :: ext
    type (status_t), intent(out) :: status

    integer :: nknots, j

    status = NF_STATUS_OK

    nknots = size(knots)

    if (x < knots(1)) then
        select case (ext)
        case (NF_INTERP_EVAL_BOUNDARY)
            ilbound = 1
            weight = 1.0_PREC
            return
        case (NF_INTERP_EVAL_CONST)
            ilbound = 0
            ! Weight does not matter in this case, should not be used
            weight = 0.0_PREC
            return
        case (NF_INTERP_EVAL_ERROR)
            status = NF_STATUS_BOUNDS_ERROR
            return
        case (NF_INTERP_EVAL_EXTRAPOLATE)
            ilbound = 1
            weight = (knots(2) - x) / (knots(2) - knots(1))
            return
        end select
    else if (x > knots(nknots)) then
        select case (ext)
        case (NF_INTERP_EVAL_BOUNDARY)
            ! Upper bound: set corresponding interpolation interval to the
            ! last one. Correct boundary value will be interpolated below.
            ilbound = nknots - 1
            weight = 0.0_PREC
            return
        case (NF_INTERP_EVAL_CONST)
            ilbound = nknots
            ! Weight does not matter in this case, should not be used
            weight = 0.0_PREC
            return
        case (NF_INTERP_EVAL_ERROR)
            status = NF_STATUS_BOUNDS_ERROR
            return
        case (NF_INTERP_EVAL_EXTRAPOLATE)
            ilbound = nknots - 1
            weight = (knots(nknots) - x) / (knots(nknots) - knots(nknots-1))
            return
        end select
    end if

    call bsearch (x, knots, j)

    ilbound = j
    weight = (knots(j+1) - x) / (knots(j+1) - knots(j))

end subroutine



pure subroutine interp_find_cached_impl_1d (x, knots, ilbound, weight, ext, &
        cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for various interpolation methods.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    integer (NF_ENUM_KIND), intent(in) :: ext
    type (search_cache), intent(inout) :: cache
    type (status_t), intent(out) :: status

    status = NF_STATUS_OK

    if (ext == NF_INTERP_EVAL_EXTRAPOLATE) then
        ! Specific implementation that assumes extrapolation
        call interp_find_impl_ext (x, knots, ilbound, weight, cache)
    else
        ! Default implementation
        call interp_find_impl_default (x, knots, ilbound, weight, ext, cache, status)
    end if

end subroutine



pure subroutine interp_find_impl_1d (x, knots, ilbound, weight, ext, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for various interpolation methods.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    integer (NF_ENUM_KIND), intent(in) :: ext
    type (status_t), intent(out) :: status

    status = NF_STATUS_OK

    if (ext == NF_INTERP_EVAL_EXTRAPOLATE) then
        ! Specific implementation that assumes extrapolation
        call interp_find_impl_ext (x, knots, ilbound, weight)
    else
        ! Default implementation
        call interp_find_impl_default (x, knots, ilbound, weight, ext, status)
    end if

end subroutine



pure subroutine interp_find_cached_impl_default_1d (x, knots, ilbound, weight, &
        ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for various interpolation methods.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    integer (NF_ENUM_KIND), intent(in) :: ext
    type (search_cache), intent(inout) :: cache
    type (status_t), intent(out) :: status

    integer :: i, nx, nknots, j
    real (PREC) :: xi, xlb, xub

    status = NF_STATUS_OK

    nknots = size(knots)
    nx = size(x)

    xlb = knots(1)
    xub = knots(nknots)

    do i = 1, nx
        xi = x(i)
        if (xi < xlb) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                ilbound(i) = 1
                weight(i) = 1.0_PREC
                cycle
            case (NF_INTERP_EVAL_CONST)
                ilbound(i) = 0
                ! Weight does not matter in this case, should not be used
                weight(i) = 0.0_PREC
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                return
            end select

            ! Default: extrapolate to the left; bracketing interval is
            ! already known, so no need to search.
            ilbound(i) = 1
            weight(i) = (knots(2) - xi) / (knots(2) - knots(1))
            cycle

        else if (xi > xub) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                ! Upper bound: set corresponding interpolation interval to the
                ! last one. Correct boundary value will be interpolated below.
                ilbound(i) = nknots - 1
                weight(i) = 0.0_PREC
                cycle
            case (NF_INTERP_EVAL_CONST)
                ilbound(i) = nknots
                ! Weight does not matter in this case, should not be used
                weight(i) = 0.0_PREC
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                return
            end select

            ! Default: extrapolate to the right; bracketing interval is
            ! already known, so no need to search.
            ilbound(i) = nknots - 1
            weight(i) = (knots(nknots) - xi) / (knots(nknots) - knots(nknots-1))
            cycle

        end if

        call bsearch_cached (xi, knots, j, cache)

        ilbound(i) = j
        weight(i) = (knots(j+1) - xi) / (knots(j+1) - knots(j))
    end do

end subroutine



pure subroutine interp_find_impl_default_1d (x, knots, ilbound, weight, ext, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for various interpolation methods.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    integer (NF_ENUM_KIND), intent(in) :: ext
    type (status_t), intent(out) :: status

    integer :: i, nx, nknots, j
    real (PREC) :: xi, xlb, xub

    status = NF_STATUS_OK

    nknots = size(knots)
    nx = size(x)

    xlb = knots(1)
    xub = knots(nknots)

    j = 1

    do i = 1, nx
        xi = x(i)
        if (xi < xlb) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                ilbound(i) = 1
                weight(i) = 1.0_PREC
                cycle
            case (NF_INTERP_EVAL_CONST)
                ilbound(i) = 0
                ! Weight does not matter in this case, should not be used
                weight(i) = 0.0_PREC
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                return
            end select

            ! Default: extrapolate to the left; bracketing interval is
            ! already known, so no need to search.
            ilbound(i) = 1
            weight(i) = (knots(2) - xi) / (knots(2) - knots(1))
            cycle

        else if (xi > xub) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                ! Upper bound: set corresponding interpolation interval to the
                ! last one. Correct boundary value will be interpolated below.
                ilbound(i) = nknots - 1
                weight(i) = 0.0_PREC
                cycle
            case (NF_INTERP_EVAL_CONST)
                ilbound(i) = nknots
                ! Weight does not matter in this case, should not be used
                weight(i) = 0.0_PREC
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                return
            end select

            ! Default: extrapolate to the right; bracketing interval is
            ! already known, so no need to search.
            ilbound(i) = nknots - 1
            weight(i) = (knots(nknots) - xi) / (knots(nknots) - knots(nknots-1))
            cycle

        end if

        call bsearch_cached_impl (xi, knots, j)

        ilbound(i) = j
        weight(i) = (knots(j+1) - xi) / (knots(j+1) - knots(j))
    end do

end subroutine



pure subroutine interp_find_cached_impl_ext_1d (x, knots, ilbound, weight, cache)
    !*  INTERP_FIND identifies the bracketing interval and interpolation
    !   weights that are required for various interpolation methods.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    type (search_cache), intent(inout) :: cache

    integer :: i, nx, nknots, j
    real (PREC) :: xi, dx, xub, wi

    nknots = size(knots)
    nx = size(x)

    ! Ensure that cache contains valid initial index as the implementation
    ! routine will not check for that!
    j = max(min(cache%i, nknots-1), 1)

    do i = 1, nx
        xi = x(i)
        call bsearch_cached_impl (xi, knots, j)

        xub = knots(j+1)
        dx = xub - knots(j)
        wi = (xub - xi) / dx

        ilbound(i) = j
        weight(i) = wi
    end do

    cache%i = j

end subroutine



pure subroutine interp_find_impl_ext_1d (x, knots, ilbound, weight)
    !*  INTERP_FIND identifies the bracketing interval and interpolation
    !   weights that are required for various interpolation methods.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight

    integer :: i, nx, j
    real (PREC) :: xi, dx, xub, wi

    nx = size(x)

    j = 1

    do i = 1, nx
        xi = x(i)
        call bsearch_cached_impl (xi, knots, j)

        xub = knots(j+1)
        dx = xub - knots(j)
        wi = (xub - xi) / dx

        ilbound(i) = j
        weight(i) = wi
    end do

end subroutine



pure subroutine interp_find_cached_1d (x, knots, ilbound, weight, ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for various interpolation methods.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    type (search_cache), intent(inout) :: cache
    type (status_t), intent(out), optional :: status

    integer (NF_ENUM_KIND) :: lext
    type (status_t) :: lstatus

    call interp_find_check_input (x, knots, ilbound, weight, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    call interp_find_impl (x, knots, ilbound, weight, lext, cache, lstatus)

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine interp_find_1d (x, knots, ilbound, weight, ext, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for various interpolation methods.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    type (status_t), intent(out), optional :: status

    integer (NF_ENUM_KIND) :: lext
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    call interp_find_check_input (x, knots, ilbound, weight, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    call interp_find_impl (x, knots, ilbound, weight, lext, lstatus)

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine interp_find_cached_0d (x, knots, ilbound, weight, ext, &
        cache, status)
    !*  INTERP_FIND_CACHE_SCALAR identifies the bracketing interval and
    !   interpolation weights that are required for various interpolation methods.
    real (PREC), intent(in) :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out) :: ilbound
    real (PREC), intent(out) :: weight
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    type (search_cache), intent(inout) :: cache
        !*  Cache that can possibly speed up search algorithm
    type (status_t), intent(out), optional :: status

    integer (NF_ENUM_KIND) :: lext
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    call interp_find_check_input (knots=knots, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    call interp_find_impl (x, knots, ilbound, weight, lext, cache, lstatus)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine interp_find_0d (x, knots, ilbound, weight, ext, status)
    !*  INTERP_FIND_SCALAR identifies the bracketing interval and
    !   interpolation weights that are required for various interpolation methods.
    real (PREC), intent(in) :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out) :: ilbound
    real (PREC), intent(out) :: weight
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    type (status_t), intent(out), optional :: status

    integer (NF_ENUM_KIND) :: lext
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    call interp_find_check_input (knots=knots, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    call interp_find_impl (x, knots, ilbound, weight, lext, status=lstatus)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine interp_find_incr_ext_impl (x, knots, ilbound, weight)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight

    integer ::nx, i, ilb, iub
    real (PREC) :: xi, wi, xub, dx, wlast

    nx = size(x)

    ! Manually compute first element
    xi = x(1)
    call bsearch (xi, knots, ilb)
    xub = knots(ilb+1)
    dx = xub - knots(ilb)
    wi = (xub - xi) / dx

    ilbound(1) = ilb
    weight(1) = wi

    if (nx == 1) return

    ! Manually compute last element
    iub = ilb
    xi = x(nx)
    call bsearch_cached_impl (xi, knots, iub)
    xub = knots(iub+1)
    dx = xub - knots(iub)
    wlast = (xub - xi) / dx

    ! Heuristically choose between binary and linear search:
    ! Do linear search if index range spans 32 or fewer elements
    ! and thus (iub - ilb + 1) = 32.
    if ((iub - ilb) <= 31) then
        do i = 2, nx - 1
            xi = x(i)
            ilb = lsearch_bounded (xi, knots, ilb, iub+1)
            xub = knots(ilb+1)
            dx = xub - knots(ilb)
            wi = (xub - xi) / dx

            ilbound(i) = ilb
            weight(i) = wi
        end do
    else
        do i = 2, nx-1
            xi = x(i)
            ilb = bsearch_bounded (xi, knots, ilb, iub+1)
            xub = knots(ilb + 1)
            dx = xub - knots(ilb)
            wi = (xub - xi) / dx

            ilbound(i) = ilb
            weight(i) = wi
        end do
    end if

    ! Write back index, weight for last element here instead of doing
    ! it immediately as elements 1 and NX could be far apart.
    ilbound(nx) = iub
    weight(nx) = wlast

end subroutine



pure subroutine interp_find_incr_ext (x, knots, ilbound, weight, status)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i

    status = NF_STATUS_OK

    call interp_find_check_input (x, knots, ilbound, weight, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Check that X is indeeded (weakly) increasing
    do i = 2, size(x)
        if (x(i) < x(i-1)) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end do

    call interp_find_incr_ext_impl (x, knots, ilbound, weight)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine interp_find_decr_ext_impl (x, knots, ilbound, weight)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight

    integer ::nx, i, ilb, iub, ilast
    real (PREC) :: xi, wi, xub, dx, wlast

    nx = size(x)

    ! Manually compute first element
    xi = x(1)
    call bsearch (xi, knots, iub)
    xub = knots(iub+1)
    dx = xub - knots(iub)
    wi = (xub - xi) / dx

    ilbound(1) = iub
    weight(1) = wi

    if (nx == 1) return

    ! Manually compute last element
    ilb = iub
    xi = x(nx)
    call bsearch_cached_impl (xi, knots, ilb)
    xub = knots(ilb+1)
    dx = xub - knots(ilb)
    wlast = (xub - xi) / dx
    ilast = ilb

    ! Heuristically choose between binary and linear search
    ! Do linear search if there are 32 elements between 1,...,NX;
    ! and thus (iub - ilb + 1) = 32.
    if ((iub - ilb) <= 31) then
        do i = 2, nx - 1
            xi = x(i)
            iub = lsearch_bounded (xi, knots, ilb, iub+1)
            xub = knots(iub+1)
            dx = xub - knots(iub)
            wi = (xub - xi) / dx

            ilbound(i) = iub
            weight(i) = wi
        end do
    else
        do i = 2, nx-1
            xi = x(i)
            iub = bsearch_bounded (xi, knots, ilb, iub+1)
            xub = knots(iub+1)
            dx = xub - knots(iub)
            wi = (xub - xi) / dx

            ilbound(i) = iub
            weight(i) = wi
        end do
    end if

    ! Write back index, weight for last element here instead of doing
    ! it immediately as elements 1 and NX could be far apart.
    ilbound(nx) = ilast
    weight(nx) = wlast

end subroutine



pure subroutine interp_find_decr_ext (x, knots, ilbound, weight, status)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i

    status = NF_STATUS_OK

    call interp_find_check_input (x, knots, ilbound, weight, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Check that X is indeeded (weakly) decreasing
    do i = 2, size(x)
        if (x(i) > x(i-1)) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end do

    call interp_find_decr_ext_impl (x, knots, ilbound, weight)

100 continue

    if (present(status)) status = lstatus

end subroutine

