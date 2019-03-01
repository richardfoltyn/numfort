


impure elemental function __APPEND(logaddexp,__PREC) (x, y) result(res)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
    real (PREC), intent(in) :: y
    real (PREC) :: res

    real (PREC) :: tmp

    real (PREC), parameter :: LOGE2 = log(2.0_PREC)

    if (x == y) then
        ! Handles infinities of the same sign without warnings
        res = x + LOGE2
    else
        tmp = x - y
        if (tmp > 0.0_PREC) then
            res = x + log1p(exp(-tmp))
        else if (tmp <= 0.0_PREC) then
            res = y + log1p(exp(tmp))
        else
            ! NaNs
            res = tmp
        end if
    end if

end function
