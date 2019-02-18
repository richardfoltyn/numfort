

pure subroutine __APPEND(ppoly_power_check_input,__PREC) (self, knots, coefs, status)
    !*  PPOLY_POWER_CHECK_INPUT implements a general-purpose input checker
    !   for piecewise polynomials using the power basis that
    !   can be used by all routines working with this data type.
    integer, parameter :: PREC = __PREC
    type (ppoly), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
    real (PREC), intent(in), dimension(:) :: coefs
    type (status_t), intent(out) :: status

    integer :: nknots, ncoefs

    nknots = ppoly_get_nknots (self)
    if (size(knots) /= nknots) goto 100

    ncoefs = ppoly_get_ncoefs (self)
    if (size(coefs) /= ncoefs) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



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
    real (PREC), intent(in), dimension(:), contiguous :: knots
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
    type (search_cache) :: bcache

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
                call bsearch_cached (xi, knots, j, bcache)
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
                call bsearch_cached (xi, knots, j, bcache)
            end select
        else
            ! Find interval j such that knots(j) <= x(i)
            call bsearch_cached (xi, knots, j, bcache)
        end if

        ! At this point we have one of three cases:
        !   1. interpolation as XLB <= X(i) <= XUB
        !   2. pseudo-interpolation that assigs the boundary value y(XLB) or y(XUB)
        !       as x(i) was outside of the interval defined by XP.
        !   3. Extrapolation for non-interior X(i) requested by user.

        ! Index of coefficient block for interval j
        jj = (j-1) * (k+1) + 1

        ! Normalize to unit interval length
        si = (xi - knots(j)) / (knots(j+1)-knots(j))
        ! Start with k-th degree term
        yi = coefs(jj+k) * si
        do ik = k-1, 1, -1
            yi = (yi + coefs(jj+ik)) * si
        end do
        ! Add constant term
        yi = yi + coefs(jj)

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
    real (PREC), intent(in), dimension(:), contiguous :: knots
        !*  x-values of data points (knots) that define the breakpoints of
        !   the piecewise polynomial
    real (PREC), intent(in), dimension(:), contiguous :: coefs
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



subroutine __APPEND(power_ppolyder,__PREC) (self, knots, coefs, poly_out, &
        coefs_out, m, status)
    !*  POWER_PPOLYDER differentiates a given polynomial and returns
    !   a new piecewise polynomial of a correpondingly lower degree.
    integer, parameter :: PREC = __PREC
    type (ppoly), intent(in) :: self
    real (PREC), intent(in), dimension(:), contiguous :: knots
        !*  Knots of the original polynomial. The differentiated polynomial
        !   will use identical knots and thus no new ones are returned.
    real (PREC), intent(in), dimension(:), contiguous, target :: coefs
        !*  Coefficient array of the original piecewise polynomial
    type (ppoly), intent(inout) :: poly_out
        !*  Differentiated polynomial
    real (PREC), intent(inout), dimension(:), contiguous, target :: coefs_out
        !*  Coefficient array of the differentiated polynomial
        !   (size must match the value returned by PPOLY_GET_NCOEFS for the
        !   given number of knots and the result polynomial degree)
    integer, intent(in), optional :: m
        !*  Order of differentiation (default: 1)
    type (status_t), intent(out), optional :: status
        !*  Optional exit code

    type (status_t) :: lstatus
    integer :: kk, deg, nknots, k1, ncoefs_out, lm, i, ncoefs, ir, ic
    real (PREC) :: dx
    real (PREC), dimension(:,:), allocatable :: tm_der
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs, ptr_coefs_out

    ! Arguments to GEMM
    integer :: m1, n, k, lda, ldb, ldc
    character (1), parameter :: transa = 'N', transb = 'N'
    real (PREC), parameter :: alpha = 1.0, beta = 0.0

    lstatus = NF_STATUS_OK

    nullify (ptr_coefs, ptr_coefs_out)

    ! Input checks
    call ppoly_power_check_input (self, knots, coefs, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_nonneg (1, m, "m", lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lm = 1
    if (present(m)) lm = m

    deg = ppoly_get_degree (self)
    nknots = ppoly_get_nknots (self)
    ncoefs = size(coefs)
    k1 = max(0, deg-lm)

    ncoefs_out = ppoly_get_ncoefs (self, k=k1)
    if (size(coefs_out) /= ncoefs_out) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    poly_out%degree = k1
    poly_out%nknots = nknots

    if (lm == 0) then
        coefs_out = coefs
        goto 100
    else if (deg == 0) then
        ! Input polynomial was a constant function, thus derivative is 0
        coefs_out = 0.0
        goto 100
    end if

    ! Build a matrix to create the coefficeints of a polynomial of deg
    ! (k-m) in a iterative faction
    k = deg
    allocate (tm_der(0:k1,0:k), source=0.0_PREC)

    if (k >= m) then
        kk = poch (k-m+1, m)
        tm_der(k-m,k) = kk
        do i = 0, k-m-1
            kk = kk / (k-i) * (k-m-i)
            ir = k - m - i - 1
            ic = k - i - 1
            tm_der(ir,ic) = kk
        end do
    end if

    ! Set up call to GEMM
    n = ncoefs/(deg+1)
    m1 = deg - lm + 1
    k = deg + 1
    lda = m1
    ldb = k
    ldc = m1

    ptr_coefs(0:deg,1:n) => coefs
    ptr_coefs_out(0:deg-m,1:n) => coefs_out

    call GEMM (transa, transb, m1, n, k, alpha, tm_der, lda, ptr_coefs, ldb, &
        beta, ptr_coefs_out, ldc)

    ! Each interval needs to be rescaled by (x_ub-x_lb)^(-m) to account
    ! for the fact that the coefficients are normalized to an interval length
    ! of 1.0
    do i = 1, n
        dx = knots(i+1) - knots(i)
        ptr_coefs_out(:,i) = ptr_coefs_out(:,i) / dx**lm
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine
