

pure subroutine __APPEND(ppolyval_check_input,__PREC) (self, knots, coefs, x, y, &
        ext, left, right, status)
    integer, parameter :: PREC = __PREC
    type (ppoly), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
    real (PREC), intent(in), dimension(:) :: coefs
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (status_t), intent(out) :: status

    integer :: ncoefs
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
    ncoefs = ppoly_get_ncoefs (self)
    if (size(coefs) /= ncoefs) then
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


pure subroutine __APPEND(ppolyval,__PREC) (self, knots, coefs, x, y, ext, &
        left, right, status)
    !*  PPOLYVAL evaluates the fitted piecewise cubic polynomial at
    !   a set of given points.
    integer, parameter :: PREC = __PREC
    type (ppoly), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
        !*  x-values of data points (knots) that define the breakpoints of
        !   the piecewise polynomial
    real (PREC), intent(in), dimension(:) :: coefs
        !*  Stacked array of polynomial coefficients for each segment
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
    integer :: i, j, jj, ik, n, nxp, k

    lstatus = NF_STATUS_OK

    call ppolyval_check_input (self, knots, coefs, x, y, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    nxp = size(knots)
    n = size(x)
    k = self%degree

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

        si = xi - knots(j)
        yi = coefs(jj)
        do ik = 1, k
            yi = yi + coefs(jj+ik) * si ** ik
        end do

        y(i) = yi

    end do

100 continue
    if (present(status)) status = lstatus

end subroutine


pure subroutine __APPEND(ppolyval_scalar,__PREC) (self, knots, coefs, x, y, &
        ext, left, right, status)
    !*  PPOLYVAL evaluates the fitted piecewise cubic polynomial at
    !   a given (scalar!) point.
    integer, parameter :: PREC = __PREC
    type (ppoly), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
        !*  x-values of data points (knots) that define the breakpoints of
        !   the piecewise polynomial
    real (PREC), intent(in), dimension(:) :: coefs
        !*  Stacked array of polynomial coefficients for each segment
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(1) :: x1, y1

    x1(1) = x
    ! Initialize to something to present unintialized variable warnings
    ! if PPOLYVAL exits with an error.
    y1 = 0.0
    call ppolyval (self, knots, coefs, x1, y1, ext, left, right, status)
    y = y1(1)

end subroutine


pure subroutine __APPEND(bernstein_ppolyval_check_input,__PREC) (self, &
        x, y, ext, left, right, status)
    !*  BERNSTEIN_PPOLYVAL_CHECK_INPUT performs input checks on caller-provided
    !   arguments to BERNSTEIN_PPOLYVAL.
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (status_t), intent(out) :: status
        !   Exit code (optional)

    if (size(x) /= size(y)) goto 100

    ! Need to provide constants for extrapolation if corresponding extrapolation
    ! mode was requested.
    if (present(ext)) then
        if (ext == NF_INTERP_EVAL_CONST .and. (.not. present(left) .or. &
            .not. present(right))) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine


pure subroutine __APPEND(bernstein_ppolyval,__PREC) (self, knots, coefs, x, y, &
        ext, left, right, status)
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
        !*  x-values of data points (knots) that define the breakpoints of
        !   the piecewise polynomial
    real (PREC), intent(in), dimension(:), contiguous :: coefs
        !*  Stacked array of polynomial coefficients for each segment
    real (PREC), intent(in), dimension(:) :: x
        !*  x-coordinates of points where function should be interpolating.
    real (PREC), intent(out), dimension(:) :: y
        !*  On exit, contains interpolated function values for x-coordinates
        !   given in X.
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Type of extrapolation, if applicable (default: extrapolate based
        !   on coefficients of closest interior interval).
    real (PREC), intent(in), optional :: left
        !*  Constant value assigned to x-coordinates smaller than KNOTS(1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    real (PREC), intent(in), optional :: right
        !*  Constant value assigned to x-coordinates larger than KNOTS(-1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    type (status_t), intent(out), optional :: status
        !   Exit code (optional)

    type (status_t) :: lstatus
    integer (NF_ENUM_KIND) :: lext

    lstatus = NF_STATUS_OK

    ! Perform common input checks
    call bernstein_check_input (self, knots, coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Perform input checks specific to evaluating polynomials
    call ppolyval_check_input (self, x, y, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    call ppolyval_impl (self, knots, coefs, x, y, lext, left, right, lstatus)

100 continue

    if (present(status)) status = lstatus

end subroutine


pure subroutine __APPEND(bernstein_ppolyval_scalar,__PREC) (self, knots, coefs, &
        x, y, ext, left, right, status)
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
        !*  x-values of data points (knots) that define the breakpoints of
        !   the piecewise polynomial
    real (PREC), intent(in), dimension(:), contiguous :: coefs
        !*  Stacked array of polynomial coefficients for each segment
    real (PREC), intent(in) :: x
        !*  x-coordinate of point where function should be interpolating.
    real (PREC), intent(out) :: y
        !*  On exit, contains interpolated function value for x-coordinate
        !   given in X.
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Type of extrapolation, if applicable (default: extrapolate based
        !   on coefficients of closest interior interval).
    real (PREC), intent(in), optional :: left
        !*  Constant value assigned to x-coordinates smaller than KNOTS(1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    real (PREC), intent(in), optional :: right
        !*  Constant value assigned to x-coordinates larger than KNOTS(-1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    type (status_t), intent(out), optional :: status
        !   Exit code (optional)

    real (PREC), dimension(1) :: x1, y1

    x1(1) = x
    ! Initialize to something to present unintialized variable warnings
    ! if POLYVAL exits with an error.
    y1 = 0.0

    call ppolyval (self, knots, coefs, x1, y1, ext, left, right, status)
    y = y1(1)

end subroutine



pure subroutine __APPEND(bernstein_ppolyval_impl,__PREC) (self, knots, coefs, &
        x, y, ext, left, right, status)
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
        !*  x-values of data points (knots) that define the breakpoints of
        !   the piecewise polynomial
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
        !*  Stacked array of polynomial coefficients for each segment
    real (PREC), intent(in), dimension(:) :: x
        !*  x-coordinates of points where function should be interpolating.
    real (PREC), intent(out), dimension(:) :: y
        !*  On exit, contains interpolated function values for x-coordinates
        !   given in X.
    integer (NF_ENUM_KIND), intent(in) :: ext
        !*  Type of extrapolation
    real (PREC), intent(in), optional :: left
        !*  Constant value assigned to x-coordinates smaller than X(1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    real (PREC), intent(in), optional :: right
        !*  Constant value assigned to x-coordinates larger than X(-1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    type (status_t), intent(out) :: status
        !   Exit code

    integer :: i, nx, k, nknots, j, jj
    real (PREC) :: xi, s, xlb, xub

    status = NF_STATUS_OK

    nknots = size(knots)
    nx = size(x)
    k = self%degree

    xlb = knots(1)
    xub = knots(nknots)

    do i = 1, nx
        xi = x(i)
        if (xi < xlb) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                j = 1
                xi = xlb
            case (NF_INTERP_EVAL_CONST)
                y(i) = left
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval j that will be used to intropolate/extrapolate
                ! y-value
                j = interp_find (xi, knots)
            end select
        else if (xi > xub) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                ! Upper bound: set corresponding interpolation interval to the
                ! last one. Correct boundary value will be interpolated below.
                j = nknots - 1
                xi = xub
            case (NF_INTERP_EVAL_CONST)
                y(i) = right
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval j such that knots(j) <= x(i)
                j = interp_find (xi, knots)
            end select
        else
            ! Find interval j such that knots(j) <= x(i)
            j = interp_find (xi, knots)
        end if

        ! Offset of coefficient block for interval j
        jj = (j-1) * (k+1)

        s = (xi - knots(j))/(knots(j+1) - knots(j))
        y(i) = bernstein_polyval (k, s, coefs(jj:jj+k))
    end do

100 continue

end subroutine


pure function __APPEND(bernstein_polyval,__PREC) (k, s, coefs) result(res)
    !*  BERNSTEIN_POLYVAL is a helper function that evaluates a polynomial
    !   wrt the Bernstein basis within a given interval.
    integer, parameter :: PREC = __PREC
    integer, intent(in) :: k
        !*  Polynomial degree
    real (PREC), intent(in) :: s
        !*  x-value relative to the interval, ie s = (x-x_lb)/(x_ub-x_lb)
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
        !*  Coefficients corresponding to given interval.
    real (PREC) :: res
        !*  Polynomial value in given interval.

    real (PREC) :: s1, comb
    integer :: j

    s1 = 1.0_PREC - s

    ! Hardcode values for lower polynomial degrees
    select case (k)
    case (0)
       res = coefs(0)
    case (1)
        res = coefs(0) * s1 + coefs(1) * s
    case (2)
        res = coefs(0) * s1**2.0_PREC + coefs(1) * 2.0_PREC*s1*s &
            + coefs(2) * s**2.0_PREC
    case (3)
        res = coefs(0) * s1**3.0 + coefs(1) * 3.0_PREC * s1**2.0_PREC * s &
            + coefs(2) * 3.0_PREC * s1 * s**2.0_PREC + coefs(3) * s**3.0_PREC
    case default
        ! TODO: Could be replaced with Casteljau's algorithm
        res = 0.0_PREC
        comb = 1.0_PREC
        do j = 0, k
            res = res + comb * s**j * s1**(k-j) * coefs(j)
            comb = comb * (k-j) / (j + 1.0_PREC)
        end do
    end select

end function