


pure subroutine __APPEND(bernstein_check_input,__PREC) (self, knots, coefs, &
        n, k, status)

    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
        !*  Array of knots
    real (PREC), intent(in), dimension(:) :: coefs
        !*  Coefficient array
    integer, intent(in), optional :: n
        !*  Number of data points. Only relevant for routines that create
        !   piecewise polynomials.
    integer, intent(in), optional :: k
        !*  Polynomial degree. Only relevant for routines that create
        !   piecewise polynomials.
    type (status_t), intent(out) :: status

    integer :: nknots, ncoefs

    nknots = ppoly_get_nknots (self, n, k)
    ncoefs = ppoly_get_ncoefs (self, n, k)

    if (size(knots) < nknots) goto 100
    if (size(coefs) /= ncoefs) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

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
        !*  x-coordinates of points where function should be interpolated
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
    type (search_cache) :: bcache

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
                call bsearch_cached (xi, knots, j, bcache)
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
                call bsearch_cached (xi, knots, j, bcache)
            end select
        else
            ! Find interval j such that knots(j) <= x(i)
            call bsearch_cached (xi, knots, j, bcache)
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

!    real (PREC), dimension(0:3) :: lcoefs3
!    integer :: i
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
!        lcoefs3(0) = coefs(0)
!        lcoefs3(1) = - 3.0_PREC * coefs(0) + 3.0_PREC * coefs(1)
!        lcoefs3(2) = 3.0_PREC * coefs(0) - 6.0_PREC * coefs(1) + 3.0_PREC * coefs(2)
!        lcoefs3(3) = -coefs(0) + 3.0_PREC * coefs(1) - 3.0_PREC * coefs(2) + coefs(3)
!        res = lcoefs3(3) * s
!        do i = 2, 1, -1
!            res = (res + lcoefs3(i)) * s
!        end do
!        res = res + lcoefs3(0)

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



pure subroutine __APPEND(bernstein_fit_deriv_check_input,__PREC) &
        (self, x, y, k, status)

    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    integer, intent(in) :: k
    type (status_t), intent(out) :: status

    integer :: n

    n = size(x)

    call check_nonneg (1, k, "k", status)
    if (status /= NF_STATUS_OK) goto 100

    ! Need at least one interval to fit a piecewise polynomial
    if (n < 2) goto 100

    if (size(x) /= size(y,2)) goto 100
    if (size(y,1) < ((k+1) / 2)) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine


pure subroutine __APPEND(bernstein_fit_deriv,__PREC) (self, x, y, k, &
        knots, coefs, status)
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(inout) :: self
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    integer, intent(in) :: k
    real (PREC), intent(out), dimension(:) :: knots
    real (PREC), intent(out), dimension(:), target, contiguous :: coefs
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: n

    lstatus = NF_STATUS_OK
    n = size(x)

    ! Perform input checking common to all routines processing Bernstein
    ! polynomials.
    call bernstein_check_input (self, knots, coefs, n, k, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Perform input checking specifing to fitting polynomials to derivatives.
    call bernstein_fit_deriv_check_input (self, x, y, k, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call bernstein_fit_deriv_impl (self, x, y, k, knots, coefs, lstatus)

100 continue
    if (present(status)) status = lstatus
end subroutine



pure subroutine __APPEND(bernstein_fit_deriv_impl,__PREC) (self, x, y, k, &
        knots, coefs, status)
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(inout) :: self
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(0:,:) :: y
    integer, intent(in) :: k
    real (PREC), intent(out), dimension(:) :: knots
    real (PREC), intent(out), dimension(:), target, contiguous :: coefs
    type (status_t), intent(out) ::  status

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs
    integer :: n, i, ka, kb, q, j, nknots
    real (PREC) :: dx, cj, p
    real (PREC), dimension(:), allocatable :: ci

    status = NF_STATUS_OK

    n = size(x)
    nknots = n
    ptr_coefs(0:k,1:n-1) => coefs

    self%degree = k
    self%nknots = nknots

    ka = size(y, 1)
    kb = ka
    allocate (ci(0:k), source=0.0_PREC)

    do i = 1, n-1
        dx = x(i+1) - x(i)
        do q = 0, ka - 1
            p = poch (k+1-q, q)
            ci(q) = y(q,i) / p * dx ** q
            cj = 0.0
            do j = 0, q-1
                cj = cj - (-1)**(j+q) * comb(q, j) * ci(j)
            end do
            ci(q) = ci(q) + cj
        end do

        do q = 0, kb - 1
            p = poch (k+1-q, q)
            ci(k-q) = y(q,i+1) / p * (-1)**q * dx**q
            cj = 0.0
            do j = 0, q-1
                cj = cj - (-1)**(j+1) * comb (q, j+1) * ci(k+1-q+j)
            end do
            ci(k-q) = ci(k-q) + cj
        end do

        ptr_coefs(:,i) = ci
    end do

    deallocate (ci)

    ! For PCHIP-type
    knots = x

end subroutine



pure subroutine __APPEND(bernstein_ppolyder,__PREC) (self, knots, coefs, &
        ppoly_out, coefs_out, m, status)
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    type (ppoly_bernstein), intent(out) :: ppoly_out
    real (PREC), intent(out), dimension(:), contiguous :: coefs_out
    integer, intent(in), optional :: m
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: k1, lm, ncoefs

    lstatus = NF_STATUS_OK

    call bernstein_check_input (self, knots, coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_nonneg (1, m, "m", lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lm = 1
    if (present(m)) lm = m

    ! TODO: Implement support for higher derivatives
    if (lm > 1) then
        lstatus = NF_STATUS_NOT_IMPLEMENTED
        goto 100
    end if

    k1 = self%degree - lm

    ncoefs = ppoly_get_ncoefs (self, k=k1)
    if (size(coefs_out) /= ncoefs) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    call ppolyder_impl (self, knots, coefs, lm, ppoly_out, coefs_out, lstatus)

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(bernstein_ppolyder_impl,__PREC) (self, knots, coefs, &
        m, ppoly_out, coefs_out, status)

    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
    integer, intent(in) :: m
    type (ppoly_bernstein), intent(out) :: ppoly_out
    real (PREC), intent(out), dimension(0:), contiguous :: coefs_out
    type (status_t), intent(out) :: status

    integer :: nknots, k1, i, j, k, ii0, ii1
    real (PREC) :: dx, dcoef

    status = NF_STATUS_OK

    nknots = self%nknots
    k = self%degree
    k1 = max(0, k - m)

    ppoly_out%degree = k1
    ppoly_out%nknots = nknots

    if (m == 0) then
        coefs_out(:) = coefs
        return
    end if

    ! Input polynomial was a constant function, thus derivative is 0
    if (k == 0) then
        coefs_out(:) = 0.0_PREC
        return
    end if

    ! Handle all *resulting* polynomial degrees >= 0.
    ii0 = 0
    ii1 = 0
    do i = 1, nknots - 1
        dx = knots(i+1) - knots(i)
        do j = 0, k1
            dcoef = coefs(ii0+j+1) - coefs(ii0+j)
            ! divide by DX to apply the chain rule as we are computing the
            ! derivative wrt. x, not s = (x-x_lb)/(x_ub-x_lb) on which the
            ! Bernstein polynomials are defined.
            coefs_out(ii1+j) = k * dcoef / dx
        end do

        ! Shift "pointer" to coefficient block corresponding to current segment
        ii0 = ii0 + (k+1)
        ii1 = ii1 + (k1+1)
    end do
end subroutine



subroutine __APPEND(bernstein2power,__PREC) (self, coefs, poly_out, coefs_out, status)
    !*  BERNSTEIN2POWER converts a piecewise polynomial with coefficients
    !   defined wrt. the Bernstein basis into a piecewise polynomial
    !   wrt. to the standard power basis.
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
        !*  Bernstein-basis polynomial
    real (PREC), intent(in), dimension(:), contiguous, target :: coefs
        !*  Coefficient array of the Bernstein-basis polynomial
    type (ppoly), intent(inout) :: poly_out
        !*  Resulting power basis polynomial
    real (PREC), intent(inout), dimension(:), contiguous, target :: coefs_out
        !*  Coefficient array of the piecewise polynomial translated onto
        !   the power basis
    type (status_t), intent(out), optional :: status
        !*  Optional exit code

    type (status_t) :: lstatus
    integer :: deg, ncoefs
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs, ptr_coefs_out

    ! Arguments to GEMM
    character (1), parameter :: transa = 'N', transb = 'N'
    real (PREC), parameter :: alpha = 1.0, beta = 0.0
    integer :: m, n, k, lda, ldb, ldc

    real (PREC), dimension(4,4), parameter :: tm3 = reshape([ integer :: &
                1, 0, 0, 0, -3, 3, 0, 0, 3, -6, 3, 0, -1, 3, -3, 1], &
            shape=[4,4], order=[2,1])
        !*  Matrix to transform polynomial deg 3 coefficients

    nullify (ptr_coefs, ptr_coefs_out)

    lstatus = NF_STATUS_OK

    deg = ppoly_get_degree (self)
    ncoefs = ppoly_get_ncoefs (self)

    if (size(coefs) /= ncoefs) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (size(coefs) /= size(coefs_out)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    poly_out%degree = self%degree
    poly_out%nknots = self%nknots

    select case (deg)
    case (0)
        ! piecewise constant function has identical coefficients for any basis
        coefs_out = coefs
    case (3)
        m = deg + 1
        n = ncoefs / (deg + 1)
        k = deg + 1
        lda = m
        ldb = k
        ldc = m

        ptr_coefs(1:k,1:n) => coefs
        ptr_coefs_out(1:m,1:n) => coefs_out

        call GEMM (transa, transb, m, n, k, alpha, tm3, lda, ptr_coefs, ldb, &
            beta, ptr_coefs_out, ldc)

    case default
        lstatus = NF_STATUS_NOT_IMPLEMENTED
        goto 100

    end select


100 continue

    if (present(status)) status = lstatus

end subroutine

