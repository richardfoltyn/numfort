


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



pure subroutine __APPEND(interp_find_impl_1d,__PREC) (x, knots, &
        ilbound, weight, ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for varios interpolation methods.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:) :: ilbound
    real (PREC), intent(out), dimension(:) :: weight
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
                goto 100
            end select
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
                goto 100
            end select
        end if

        call bsearch_cached (xi, knots, j, cache)

        ilbound(i) = j
        weight(i) = (knots(j+1) - xi) / (knots(j+1) - knots(j))
    end do

100 continue

end subroutine



pure subroutine __APPEND(interp_find_1d,__PREC) (x, knots, &
        ilbound, weight, ext, cache, status)
    !*  INTERP_LOCATE identifies the bracketing interval and interpolation
    !   weights that are required for varios interpolation methods.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(out), dimension(:) :: ilbound
    real (PREC), intent(out), dimension(:) :: weight
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

    real (PREC), dimension(1) :: x1d, weight1d
    integer, dimension(1) :: ilbound1d

    x1d(1) = x

    call interp_find_impl (x1d, knots, ilbound1d, weight1d, ext, cache, &
        status)

    if (status == NF_STATUS_OK) then
        ilbound = ilbound1d(1)
        weight = weight1d(1)
    end if

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
