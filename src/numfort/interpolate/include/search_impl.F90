


pure function __APPEND(bsearch,__PREC) (needle, haystack) result(lb)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: needle
    real (PREC), intent(in), dimension(:) :: haystack

    integer :: lb, ub, midx

    lb = 1
    ub = size(haystack)

    do while (ub > (lb + 1))
        midx = (ub + lb) / 2
        if (haystack(midx) > needle) then
            ub = midx
        else
            lb = midx
        end if
    end do

end function


pure function __APPEND(interp_find,__PREC) (needle, haystack) result (res)
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32
    real (PREC), intent(in), dimension(:) :: haystack
    real (PREC), intent(in) :: needle

    integer (INTSIZE) :: res, n

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

