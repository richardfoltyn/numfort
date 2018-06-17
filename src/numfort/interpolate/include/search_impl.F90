


pure function __APPEND(bsearch,__PREC) (needle, haystack, i) result(ilb)
    !*  BSEARCH finds the interval on array HAYSTACK that contains NEEDLE,
    !   returning the array index of the lower boundary of that interval.
    !
    !   If NEEDLE is smaller than the first element of HAYSTACK, the rountine
    !   returns 1. If NEEDLE is larger than the last of element of HAYSTACK,
    !   size(HAYSTACK)-1 is returned.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:) :: haystack
    integer, intent(in), optional :: i

    integer :: ilb, iub, imid, n, ii

    n = size(haystack)

    ilb = min(1, n)
    if (n <= 1) return

    if (present(i)) then
        ii = i
        ii = max(min(ii, n), 1)
        if (haystack(ii) <= needle) then
            ilb = ii
            if (ii == (n-1)) then
                return
            else if (haystack(ii+1) > needle) then
                return
            end if
            iub = n
        else
            ilb = 1
            iub = ii
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

end function


pure function __APPEND(interp_find,__PREC) (needle, haystack, i) result (res)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: haystack
    real (PREC), intent(in) :: needle
    integer, intent(in), optional :: i

    integer :: res, n

    n = size(haystack)

    if (n < 2) then
        res = -1
    else if (needle <= haystack(1)) then
        res = 1
    else if (needle >= haystack(n-1)) then
        res = n-1
    else
        res = bsearch (needle, haystack, i)
    end if

end function

