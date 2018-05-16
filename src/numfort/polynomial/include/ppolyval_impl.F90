

pure subroutine __APPEND(ppolyval_input_check,__PREC) (knots, coef, k, x, y, &
        ext, left, right, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: knots
    real (PREC), intent(in), dimension(:) :: coef
    integer, intent(in) :: k
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (status_t), intent(out) :: status

    integer :: n
    status = NF_STATUS_OK

    if (size(x) /= size(y)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Check min. data point array size
    if (size(knots) < 2) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Check whether coefficient array is at least as large as expected to
    ! prevent possible out-of-bounds access
    n = ppoly_get_ncoef (size(knots), k)
    if (size(coef) < n) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Need to provide constants for extrapolation if corresponding extrapolation
    ! mode was requested.
    if (present(ext)) then
        if (ext == NF_INTERP_EVAL_CONST .and. (.not. present(left) .or. &
            .not. present(right))) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

100 continue

end subroutine


pure subroutine __APPEND(ppolyval,__PREC) (knots, coef, k, x, y, ext, &
        left, right, status)
    !*  PPOLYVAL evaluates the fitted piecewise cubic polynomial at
    !   a set of given points.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: knots
        !*  x-values of data points (knots) that define the breakpoints of
        !   the piecewise polynomial
    real (PREC), intent(in), dimension(:) :: coef
        !*  Stacked array of polynomial coefficients for each segment
    integer, intent(in) :: k
        !*  Polynomial degree on each segment
    real (PREC), intent(in), dimension(:) :: x
        !*  x-coordinates of points where function should be interpolating.
    real (PREC), intent(out), dimension(:) :: y
        !*  On exit, contains interpolated function values for x-coordinates
        !   given in X.
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Type of extrapolation, if applicable (default: extrapolate based
        !   on coefficients of closest interior interval).
    real (PREC), intent(in), optional :: left
        !*  Constant value assigned to x-coordinates smaller than X(1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    real (PREC), intent(in), optional :: right
        !*  Constant value assigned to x-coordinates larger than X(-1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    type (status_t), intent(out), optional :: status
        !   Exit code (optional)

    type (status_t) :: lstatus
    integer (NF_ENUM_KIND) :: lext
    real (PREC) :: xlb, xub, xi, yi, si
    integer :: i, j, jj, ik, n, nxp

    lstatus = NF_STATUS_OK

    call ppolyval_input_check (knots, coef, k, x, y, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    nxp = size(knots)
    n = size(x)

    xlb = knots(1)
    xub = knots(nxp)

    do i = 1, n
        xi = x(i)
        if (xi < xlb) then
            select case (lext)
            case (NF_INTERP_EVAL_BOUNDARY)
                ! Lower bound: set corresponding interpolation interval to the
                ! first one. Correct boundary value will be interpolated below.
                j = 1
                xi = xlb
            case (NF_INTERP_EVAL_CONST)
                y(i) = left
                cycle
            case (NF_INTERP_EVAL_ERROR)
                lstatus = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval j such that knots(j) <= x(i)
                j = interp_find (x(i), knots)
            end select
        else if (xi > xub) then
            select case (lext)
            case (NF_INTERP_EVAL_BOUNDARY)
                ! Upper bound: set corresponding interpolation interval to the
                ! last one. Correct boundary value will be interpolated below.
                j = nxp - 1
                xi = xub
            case (NF_INTERP_EVAL_CONST)
                y(i) = right
                cycle
            case (NF_INTERP_EVAL_ERROR)
                lstatus = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval j such that knots(j) <= x(i)
                j = interp_find (x(i), knots)
            end select
        else
            ! Find interval j such that knots(j) <= x(i)
            j = interp_find (x(i), knots)
        end if

        ! At this point we have one of three cases:
        !   1. interpolation as XLB <= X(i) <= XUB
        !   2. pseudo-interpolation that assigs the boundary value y(XLB) or y(XUB)
        !       as x(i) was outside of the interval defined by XP.
        !   3. Extrapolation for non-interior X(i) requested by user.

        ! Index of coefficient block for interval j
        jj = (j-1) * (k+1) + 1

        yi = 0.0
        si = xi - knots(j)
        do ik = 0, k
            yi = yi +  coef(jj+ik) * si ** ik
        end do

        y(i) = yi

    end do

100 continue
    if (present(status)) status = lstatus

end subroutine


pure subroutine __APPEND(ppolyval_scalar,__PREC) (knots, coef, k, x, y, &
        ext, left, right, status)
    !*  PPOLYVAL evaluates the fitted piecewise cubic polynomial at
    !   a given (scalar!) point.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: knots
        !*  x-values of data points (knots) that define the breakpoints of
        !   the piecewise polynomial
    real (PREC), intent(in), dimension(:) :: coef
        !*  Stacked array of polynomial coefficients for each segment
    integer, intent(in) :: k
        !*  Polynomial degree on each segment
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(1) :: x1, y1
    type (status_t) :: lstatus

    x1(1) = x
    call ppolyval (knots, coef, k, x1, y1, ext, left, right, lstatus)

    ! Note: y1 may be uninitialized if error is encountered and evaluation is
    ! is not performed.
    if (lstatus == NF_STATUS_OK) then
        y = y1(1)
    end if

    if (present(status)) status = lstatus

end subroutine
