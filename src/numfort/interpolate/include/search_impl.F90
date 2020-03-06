


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



pure subroutine __APPEND(interp_find_check_input,__PREC) (x, knots, ilbound, &
        weight, status)
    !*  INTERP_FIND_CHECK_INPUT performs input validation for various
    !   INTERP_FIND routines.

    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), optional :: x
        !*  Array of points for which to determine the location within `knots`
    real (PREC), intent(in), dimension(:) :: knots
        !*  Array defining the points over which to search.
    integer (INTSIZE), intent(in), dimension(:), optional :: ilbound
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



pure subroutine __APPEND(interp_find_impl_scalar,__PREC) (x, knots, &
        ilbound, weight, ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for varios interpolation methods.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out) :: ilbound
    real (PREC), intent(out) :: weight
    integer (NF_ENUM_KIND), intent(in) :: ext
    type (search_cache), intent(inout) :: cache
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



pure subroutine __APPEND(interp_find_impl_1d,__PREC) (x, knots, &
        ilbound, weight, ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for varios interpolation methods.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    integer (NF_ENUM_KIND), intent(in) :: ext
    type (search_cache), intent(inout) :: cache
    type (status_t), intent(out) :: status

    if (ext == NF_INTERP_EVAL_EXTRAPOLATE) then
        call interp_find_impl_ext (x, knots, ilbound, weight, ext, cache, status)
    else
        call interp_find_impl_default (x, knots, ilbound, weight, ext, cache, status)
    end if

end subroutine



pure subroutine __APPEND(interp_find_impl_default_1d,__PREC) (x, knots, &
        ilbound, weight, ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for varios interpolation methods.
    integer, parameter :: PREC = __PREC
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



pure subroutine __APPEND(interp_find_impl_ext_1d,__PREC) (x, knots, &
        ilbound, weight, ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for varios interpolation methods.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    integer (NF_ENUM_KIND), intent(in) :: ext
    type (search_cache), intent(inout) :: cache
    type (status_t), intent(out) :: status

    integer :: i, nx, nknots, j
    real (PREC) :: xi, dx, const, slope

    status = NF_STATUS_OK

    nknots = size(knots)
    nx = size(x)

    do i = 1, nx
        call bsearch_cached (x(i), knots, ilbound(i), cache)
    end do

    do i = 1, nx
        xi = x(i)
        j = ilbound(i)
        dx = knots(j+1) - knots(j)
        const = knots(j+1) / dx
        slope = - 1.0_PREC / dx
        weight(i) = slope * xi + const
    end do
end subroutine



pure subroutine __APPEND(interp_find_1d,__PREC) (x, knots, &
        ilbound, weight, ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for varios interpolation methods.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:), contiguous :: ilbound
    real (PREC), intent(out), dimension(:), contiguous :: weight
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    type (search_cache), intent(inout), optional :: cache
    type (status_t), intent(out), optional :: status

    type (search_cache) :: lcache
    integer (NF_ENUM_KIND) :: lext
    type (status_t) :: lstatus

    call interp_find_check_input (x, knots, ilbound, weight, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext
    if (present(cache)) lcache = cache

    call interp_find_impl (x, knots, ilbound, weight, lext, lcache, lstatus)

    if (present(cache)) cache = lcache

100 continue
    if (present(status)) status = lstatus

end subroutine





pure subroutine __APPEND(interp_find_scalar,__PREC) (x, knots, &
        ilbound, weight, ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for varios interpolation methods.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out) :: ilbound
    real (PREC), intent(out) :: weight
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    type (search_cache), intent(inout), optional :: cache
    type (status_t), intent(out), optional :: status

    type (search_cache) :: lcache
    integer (NF_ENUM_KIND) :: lext
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    call interp_find_check_input (knots=knots, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext
    if (present(cache)) lcache = cache

    call interp_find_impl (x, knots, ilbound, weight, lext, lcache, lstatus)

    if (present(cache)) cache = lcache

100 continue

    if (present(status)) status = lstatus

end subroutine
