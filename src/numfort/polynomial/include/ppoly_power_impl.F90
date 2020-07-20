

pure subroutine __APPEND(power_ppoly_check_input,__PREC) (self, knots, coefs, &
        ilbound, weight, x, y, status)
    !*  POWER_PPOLY_CHECK_INPUT implements a general-purpose input checker
    !   for piecewise polynomials using the power basis that
    !   can be used by all routines working with this data type.
    integer, parameter :: PREC = __PREC
    type (ppoly), intent(in) :: self
    real (PREC), intent(in), dimension(:), optional :: knots
        !*  Optional array of knots. Array not required for EVAL routines.
    real (PREC), intent(in), dimension(:) :: coefs
        !*  Array of polynomial coefficients.
    integer, intent(in), dimension(:), optional :: ilbound
        !*  Array of indices of bracket lower bounds. Only required
        !   by EVAL routines.
    real (PREC), intent(in), dimension(:), optional :: weight
        !*  Interpolating weight for each bracket lower bound. Only required
        !   by EVAL functions.
    real (PREC), intent(in), dimension(:), optional :: x, y
        !*  Input/output arrays containing the points at which polynomial
        !   should be evaluated, and used to store polynomial values.
        !   Only required for 1d-array API.
    type (status_t), intent(out) :: status

    integer :: nknots, ncoefs

    if (present(knots)) then
        ! Need at least one bracket
        if (size(knots) < 2) goto 100

        nknots = ppoly_get_nknots (self)
        if (size(knots) /= nknots) goto 100
    end if

    ncoefs = ppoly_get_ncoefs (self)
    if (size(coefs) /= ncoefs) goto 100

    if (present(ilbound) .and. present(weight)) then
        if (size(ilbound) /= size(weight)) goto 100
    end if

    if (present(x) .and. present(y)) then
        if (size(x) /= size(y)) goto 100
    end if

    if (present(ilbound) .and. present(y)) then
        if (size(ilbound) /= size(y)) goto 100
    end if

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



pure subroutine __APPEND(power_ppolyval_eval_impl_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)
    !*  PPOLYVAL evaluates a piecewise polynomial at a set of given points.
    !
    !   This routine contains only the core implementation, no input validation
    !   is performed.

    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly), intent(in) :: self
    integer (INTSIZE), intent(in), dimension(:), contiguous :: ilbound
    real (PREC), intent(in), dimension(:), contiguous :: weight
    real (PREC), intent(in), dimension(:), contiguous :: coefs
        !*  Stacked array of polynomial coefficients for each segment
    real (PREC), intent(out), dimension(:), contiguous :: y
        !*  On exit, contains interpolated function values for x-coordinates
        !   given in X.
    integer (NF_ENUM_KIND), intent(in) :: ext
        !*  Type of extrapolation, if applicable (default: extrapolate based
        !   on coefficients of closest interior interval).
    real (PREC), intent(in), optional :: left
        !*  Constant value assigned to x-coordinates smaller than X(1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    real (PREC), intent(in), optional :: right
        !*  Constant value assigned to x-coordinates larger than X(-1), if
        !   EXT = NF_INTERP_EVAL_CONST.

    select case (self%degree)
    case (2)
        call ppolyval_eval_impl_deg2 (self, ilbound, weight, coefs, y, ext, left, right)
    case (3)
        call ppolyval_eval_impl_deg3 (self, ilbound, weight, coefs, y, ext, left, right)
    case default
        call ppolyval_eval_impl_degk (self, ilbound, weight, coefs, y, ext, left, right)
    end select

end subroutine



pure subroutine __APPEND(power_ppolyval_eval_impl_scalar,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)

    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC
    type (ppoly), intent(in) :: self
    integer (INTSIZE), intent(in) :: ilbound
    real (PREC), intent(in) :: weight
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
    !*  Stacked array of polynomial coefficients for each segment
    real (PREC), intent(out) :: y
    !*  On exit, contains interpolated function value for x-coordinate
    !   given in X.
    integer (NF_ENUM_KIND), intent(in) :: ext
    !*  Type of extrapolation, if applicable (default: extrapolate based
    !   on coefficients of closest interior interval).
    real (PREC), intent(in), optional :: left
    !*  Constant value assigned to x-coordinates smaller than KNOTS(1), if
    !   EXT = NF_INTERP_EVAL_CONST.
    real (PREC), intent(in), optional :: right
    !*  Constant value assigned to x-coordinates larger than KNOTS(-1), if
    !   EXT = NF_INTERP_EVAL_CONST.

    integer :: k, jj
    real (PREC) :: s

    if (ext == NF_INTERP_EVAL_CONST) then
        if (ilbound == 0) then
            y = left
            return
        else if (ilbound == self%nknots) then
            y = right
            return
        end if
    end if

    ! Evaluation point rescaled onto [0, 1]
    s = 1.0_PREC - weight

    k = self%degree
    jj = (ilbound - 1) * (k+1)
    y = power_polyval (self%degree, s, coefs(jj:jj+k))

end subroutine



pure subroutine __APPEND(power_ppolyval_eval_impl_degk_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)
    !*  POWER_PPOLYVAL_EVAL_IMPL_DEGK evaluates a piecewise polynomial wrt.
    !   the power basis of degree K at a set of given points.
    !
    !   (Internal) implementation routine, does not perform any input validation.
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly), intent(in) :: self
    integer, intent(in), dimension(:), contiguous :: ilbound
    real (PREC), intent(in), dimension(:), contiguous :: weight
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer (NF_ENUM_KIND), intent(in) :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right

    integer (INTSIZE) :: i, nx, nknots, jj, k, ilb
    real (PREC) :: s, yi, wgt

    k = self%degree
    nx = size(y)

    if (ext == NF_INTERP_EVAL_CONST) then
        ! Deal with constant "extrapolation" separately
        nknots = self%nknots

        do i = 1, nx
            ilb = ilbound(i)
            if (ilb == 0) then
                y(i) = left
            else if (ilb == nknots) then
                y(i) = right
            else
                wgt = weight(i)
                ! For interior points, position of x within the normalized
                ! unit interval is 1 - weight on bracket lower bound
                s = 1.0_PREC - wgt

                ! Offset of coefficient block for interval j
                jj = (ilb-1) * (k+1)

                yi = power_polyval (k, s, coefs(jj:jj+k))

                y(i) = yi
            end if
        end do
    else
        ! Extrapolation taken care of by adjusted weights  and bracket lower
        ! bound indices, no further adjustment required.
        do i = 1, nx
            ilb = ilbound(i)
            wgt = weight(i)

            ! For interior points, position of x within the normalized
            ! unit interval is 1 - weight on bracket lower bound
            s = 1.0_PREC - wgt

            ! Offset of coefficient block for interval j
            jj = (ilb-1) * (k+1)

            yi = power_polyval (k, s, coefs(jj:jj+k))

            y(i) = yi
        end do
    end if

end subroutine



pure subroutine __APPEND(power_ppolyval_eval_impl_deg2_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)
    !*  POWER_PPOLYVAL_EVAL_IMPL_DEG2 evaluates a piecewise quadratic polynomials
    !   wrt. the power basis at a set of given points.
    !
    !   (Internal) implementation routine, does not perform any input validation.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32
    type (ppoly), intent(in) :: self
    integer (INTSIZE), intent(in), dimension(:), contiguous :: ilbound
    real (PREC), intent(in), dimension(:), contiguous :: weight
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer (NF_ENUM_KIND), intent(in) :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right

    integer, parameter :: k = 2
    integer (INTSIZE) :: i, nx, nknots, jj, ilb
    real (PREC) :: s, yi, wgt

    nx = size(y)

    if (ext == NF_INTERP_EVAL_CONST) then
        ! Deal with constant "extrapolation" separately
        nknots = self%nknots

        do i = 1, nx
            ilb = ilbound(i)
            if (ilb == 0) then
                y(i) = left
            else if (ilb == nknots) then
                y(i) = right
            else
                wgt = weight(i)
                ! For interior points, position of x within the normalized
                ! unit interval is 1 - weight on bracket lower bound
                s = 1.0_PREC - wgt

                ! Offset of coefficient block for interval j
                jj = (ilb-1) * (k+1)

                yi = power_polyval (k, s, coefs(jj:jj+k))

                y(i) = yi
            end if
        end do
    else
        ! Extrapolation taken care of by adjusted weights  and bracket lower
        ! bound indices, no further adjustment required.
        do i = 1, nx
            ilb = ilbound(i)
            wgt = weight(i)

            ! For interior points, position of x within the normalized
            ! unit interval is 1 - weight on bracket lower bound
            s = 1.0_PREC - wgt

            ! Offset of coefficient block for interval j
            jj = (ilb-1) * (k+1)

            yi = coefs(jj+2) * s + coefs(jj+1)
            yi = yi * s + coefs(jj)

            y(i) = yi
        end do
    end if

end subroutine



pure subroutine __APPEND(power_ppolyval_eval_impl_deg3_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)
    !*  POWER_PPOLYVAL_EVAL_IMPL_DEG3 evaluates a piecewise cubic polynomials
    !   wrt. the power basis at a set of given points.
    !
    !   (Internal) implementation routine, does not perform any input validation.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32
    type (ppoly), intent(in) :: self
    integer (INTSIZE), intent(in), dimension(:), contiguous :: ilbound
    real (PREC), intent(in), dimension(:), contiguous :: weight
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer (NF_ENUM_KIND), intent(in) :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right

    integer, parameter :: k = 3
    integer (INTSIZE) :: i, j, nx, nknots, jj, ilb
    real (PREC) :: s, yi, wgt

    nx = size(y)

    if (ext == NF_INTERP_EVAL_CONST) then
        ! Deal with constant "extrapolation" separately
        nknots = self%nknots

        do i = 1, nx
            ilb = ilbound(i)
            if (ilb == 0) then
                y(i) = left
            else if (ilb == nknots) then
                y(i) = right
            else
                wgt = weight(i)
                ! For interior points, position of x within the normalized
                ! unit interval is 1 - weight on bracket lower bound
                s = 1.0_PREC - wgt

                ! Offset of coefficient block for interval j
                jj = (ilb-1) * (k+1)

                yi = power_polyval (k, s, coefs(jj:jj+k))

                y(i) = yi
            end if
        end do
    else
        ! Extrapolation taken care of by adjusted weights  and bracket lower
        ! bound indices, no further adjustment required.
        do i = 1, nx
            ilb = ilbound(i)
            wgt = weight(i)

            ! For interior points, position of x within the normalized
            ! unit interval is 1 - weight on bracket lower bound
            s = 1.0_PREC - wgt

            ! Offset of coefficient block for interval j
            jj = (ilb-1) * (k+1)

            yi = coefs(jj+3) * s + coefs(jj+2)
            do j = 1, 0, -1
                yi = yi * s + coefs(jj+j)
            end do

            y(i) = yi
        end do
    end if

end subroutine



pure subroutine __APPEND(power_ppolyval_eval_scalar,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right, status)
    !*  POWER_PPOLYVAL_EVAL implements the user-friendly front-end to
    !   the evaluation of polynomials using the power basis.
    !
    !   This routine performs input validation on given arguments.
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly), intent(in) :: self
    integer (INTSIZE), intent(in) :: ilbound
        !*  Index of bracketing interval lower bound
    real (PREC), intent(in) :: weight
        !*  Weight on bracketing interval lower bound
    real (PREC), intent(in), dimension(:), contiguous :: coefs
        !*  Piecewise polynomial coefficients
    real (PREC), intent(out) :: y
        !*  Polynomial value at given point
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer (NF_ENUM_KIND) :: lext

    lstatus = NF_STATUS_OK

    call check_input_ext (1.0_PREC, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    call power_ppoly_check_input (self, coefs=coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Evaluate polynomial at given location
    call ppolyval_eval_impl (self, ilbound, weight, coefs, y, lext, left, right)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(power_ppolyval_eval_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right, status)
    !*  POWER_PPOLYVAL_EVAL implements the user-friendly front-end to
    !   the evaluation of polynomials using the power basis.
    !   This routine performs input validation on given arguments.
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly), intent(in) :: self
    integer (INTSIZE), intent(in), dimension(:), contiguous :: ilbound
    !*  Array of indices of bracketing interval lower bounds, one for
    !   each point at which polynomial should be evaluated.
    real (PREC), intent(in), dimension(:), contiguous :: weight
    !*  Array of weights on bracketing interval lower bounds, one for
    !   each point at which polynomial should be evaluated.
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    !*  Array of piecewise polynomial coefficients
    real (PREC), intent(out), dimension(:), contiguous :: y
    !*  Polynomial values at given points.
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer (NF_ENUM_KIND) :: lext

    lstatus = NF_STATUS_OK

    call check_input_ext (1.0_PREC, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    call power_ppoly_check_input (self, coefs=coefs, ilbound=ilbound, &
        weight=weight, y=y, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Evaluate polynomial at given location
    call ppolyval_eval_impl (self, ilbound, weight, coefs, y, lext, left, right)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(power_ppolyval_scalar,__PREC) &
        (self, knots, coefs, x, y, ext, left, right, cache, status)
    !*  POWER_PPOLYVAL_SCALAR computes a piecewise polynomial defined
    !   wrt the power basis at a given point.
    !
    !   User-friedly routine that performs input checks, locates the relevant
    !   bracketing interval and evaluates the polynomial.
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly), intent(in) :: self
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (search_cache), intent(inout), optional :: cache
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (search_cache) :: lcache
    integer (NF_ENUM_KIND) :: lext
    integer (INTSIZE) :: ilbound
    real (PREC) :: weight

    lstatus = NF_STATUS_OK

    call check_input_ext (1.0_PREC, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    call power_ppoly_check_input (self, knots, coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(cache)) lcache = cache

    ! Find bracket
    call interp_find_impl (x, knots, ilbound, weight, lext, lcache, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Evaluate polynomial at given location
    call ppolyval_eval_impl (self, ilbound, weight, coefs, y, lext, left, right)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(power_ppolyval_1d,__PREC) &
        (self, knots, coefs, x, y, ext, left, right, cache, status)
    !*  POWER_PPOLYVAL computes a piecewise polynomial defined
    !   wrt the power basis for a given set of point.
    !
    !   User-friedly routine that performs input checks, locates the relevant
    !   bracketing interval for each point and evaluates the polynomial.
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly), intent(in) :: self
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (search_cache), intent(inout), optional :: cache
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (search_cache) :: lcache
    integer (NF_ENUM_KIND) :: lext
    integer (INTSIZE), dimension(:), allocatable :: ilbound
    real (PREC), dimension(:), allocatable :: weight
    integer :: nx

    lstatus = NF_STATUS_OK

    call check_input_ext (1.0_PREC, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    call power_ppoly_check_input (self, knots, coefs, x=x, y=y, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(cache)) lcache = cache

    nx = size(x)
    allocate (ilbound(nx), weight(nx))

    ! Find bracket
    call interp_find_impl (x, knots, ilbound, weight, lext, lcache, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Evaluate polynomial at given location
    call ppolyval_eval_impl (self, ilbound, weight, coefs, y, lext, left, right)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure function __APPEND(power_polyval,__PREC) (k, s, coefs) result(res)
    !*  POWER_POLYVAL evaluates a polynomial wrt the power basis at a single
    !   given point.
    integer, parameter :: PREC = __PREC

    integer, intent(in) :: k
        !*  Polynomial degree
    real (PREC), intent(in) :: s
        !*  Point at which polynomial should be evaluated, expressed relative
        !   to the bracketing interval as a value in [0, 1]
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
        !*  Polynomial coefficients for given bracketing interval, arranged
        !   in increasing exponent order.
    real (PREC) :: res
        !*  Polynomial value at given point.

    integer :: j

    res = coefs(k)
    do j = k-1, 0, -1
        res = res * s + coefs(j)
    end do

end function



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
    call power_ppoly_check_input (self, knots, coefs, status=lstatus)
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

    if (k >= lm) then
        kk = poch (k-lm+1, lm)
        tm_der(k-lm,k) = kk
        do i = 0, k-lm-1
            kk = kk / (k-i) * (k-lm-i)
            ir = k - lm - i - 1
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
    ptr_coefs_out(0:deg-lm,1:n) => coefs_out

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
