


pure subroutine __APPEND(bernstein_check_input,__PREC) (self, knots, coefs, &
        n, k, status)

    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:), optional :: knots
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

    if (present(knots)) then
        ! Knots are optional, need not be passed to EVAL routines
        if (size(knots) < nknots) goto 100
    end if

    if (size(coefs) /= ncoefs) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



pure subroutine __APPEND(bernstein_ppolyval_eval_impl_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
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
        !*  Constant value assigned to x-coordinates smaller than KNOTS(1), if
        !   EXT = NF_INTERP_EVAL_CONST.
    real (PREC), intent(in), optional :: right
        !*  Constant value assigned to x-coordinates larger than KNOTS(-1), if
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



pure subroutine __APPEND(bernstein_ppolyval_eval_impl_scalar,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    integer (INTSIZE), intent(in) :: ilbound
    real (PREC), intent(in) :: weight
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
        !*  Stacked array of polynomial coefficients for each segment
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

    s = 1.0_PREC - weight

    k = self%degree
    jj = (ilbound - 1) * (k+1)
    y = bernstein_polyval (self%degree, s, coefs(jj:jj+k))

end subroutine



pure subroutine __APPEND(bernstein_ppolyval_eval_impl_degk_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)
    !*  BERNSTEIN_PPOLYVAL_EVAL_IMPL_DEGN implements the evaluation of
    !   piecewise polynomials or arbitrary degree wrt. the Bernstein basis.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32
    type (ppoly_bernstein), intent(in) :: self
    integer (INTSIZE), intent(in), dimension(:), contiguous :: ilbound
    real (PREC), intent(in), dimension(:), contiguous :: weight
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
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

                yi = bernstein_polyval (k, s, coefs(jj:jj+k))

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

            yi = bernstein_polyval (k, s, coefs(jj:jj+k))

            y(i) = yi
        end do
    end if

end subroutine




pure subroutine __APPEND(bernstein_ppolyval_eval_impl_deg2_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)
    !*  BERNSTEIN_PPOLYVAL_EVAL_IMPL_DEG2 implements the evaluation of
    !   piecewise quadratic polynomials wrt. the Bernstein basis.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32
    type (ppoly_bernstein), intent(in) :: self
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

                yi = bernstein_polyval (k, s, coefs(jj:jj+k))

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

            yi = bernstein_polyval (k, s, coefs(jj:jj+k))

            y(i) = yi
        end do
    end if

end subroutine



pure subroutine __APPEND(bernstein_ppolyval_eval_impl_deg3_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right)
    !*  BERNSTEIN_PPOLYVAL_IMPL_CUBIC implements the evaluation of
    !   piecewise cubic polynomials wrt. the Bernstein basis.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32
    type (ppoly_bernstein), intent(in) :: self
    integer (INTSIZE), intent(in), dimension(:), contiguous :: ilbound
    real (PREC), intent(in), dimension(:), contiguous :: weight
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer (NF_ENUM_KIND), intent(in) :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right

    integer, parameter :: k = 3
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

                yi = bernstein_polyval (k, s, coefs(jj:jj+k))

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

            yi = bernstein_polyval (k, s, coefs(jj:jj+k))

            y(i) = yi
        end do
    end if

end subroutine



pure subroutine __APPEND(bernstein_ppolyval_eval_scalar,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right, status)
    !*  BERNSTEIN_PPOLYVAL_EVAL implements the user-friendly front-end to
    !   the evaluation of polynomials using the Bernstein basis.
    !   This routine performs input validation on given arguments.
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly_bernstein), intent(in) :: self
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

    call bernstein_check_input (self, coefs=coefs, n=self%nknots, k=self%degree, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Evaluate polynomial at given location
    call ppolyval_eval_impl (self, ilbound, weight, coefs, y, lext, left, right)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(bernstein_ppolyval_eval_1d,__PREC) &
        (self, ilbound, weight, coefs, y, ext, left, right, status)
    !*  BERNSTEIN_PPOLYVAL_EVAL implements the user-friendly front-end to
    !   the evaluation of polynomials using the Bernstein basis.
    !   This routine performs input validation on given arguments.
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly_bernstein), intent(in) :: self
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

    call bernstein_check_input (self, coefs=coefs, n=self%nknots, k=self%degree, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (size(ilbound) /= size(weight) .or. size(ilbound) /= size(y)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Evaluate polynomial at given location
    call ppolyval_eval_impl (self, ilbound, weight, coefs, y, lext, left, right)

    100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(bernstein_ppolyval_scalar,__PREC) &
        (self, knots, coefs, x, y, ext, left, right, cache, status)

    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly_bernstein), intent(in) :: self
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

    call bernstein_check_input (self, knots, coefs, self%nknots, self%degree, lstatus)
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



pure subroutine __APPEND(bernstein_ppolyval_1d,__PREC) &
        (self, knots, coefs, x, y, ext, left, right, cache, status)

    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    type (ppoly_bernstein), intent(in) :: self
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

    call bernstein_check_input (self, knots, coefs, self%nknots, self%degree, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Additionally check for conformable array sizes
    if (size(x) /= size(y)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

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

    real (PREC), dimension(0:3) :: lcoefs3
    real (PREC) :: s1, comb
    integer :: j

    ! Hardcode values for lower polynomial degrees
    select case (k)
    case (0)
       res = coefs(0)
    case (1)
        s1 = 1.0_PREC - s
        res = coefs(0) * s1 + coefs(1) * s

    case (2)
        lcoefs3(0) = coefs(0)
        lcoefs3(1) = 2.0_PREC * (coefs(1) - coefs(0))
        lcoefs3(2) = coefs(0) - 2.0_PREC * coefs(1) + coefs(2)

        res = lcoefs3(2) * s + lcoefs3(1)
        res = res * s + lcoefs3(0)
    case (3)
        lcoefs3(0) = coefs(0)
        lcoefs3(1) = - 3.0_PREC * coefs(0) + 3.0_PREC * coefs(1)
        lcoefs3(2) = 3.0_PREC * coefs(0) - 6.0_PREC * coefs(1) + 3.0_PREC * coefs(2)
        lcoefs3(3) = -coefs(0) + 3.0_PREC * coefs(1) - 3.0_PREC * coefs(2) + coefs(3)

        res = lcoefs3(3) * s + lcoefs3(2)
        do j = 1, 0, -1
            res = res * s + lcoefs3(j)
        end do

    case default
        ! TODO: Could be replaced with Casteljau's algorithm
        s1 = 1.0_PREC - s
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


!
!pure subroutine __APPEND(bernstein2power_deg3,__PREC) (coefs_bb, coefs_pb)
!
!    integer, parameter :: PREC = __PREC
!    real (PREC), intent(in), dimension(0:), contiguous :: coefs_bb
!    real (PREC), intent(out), dimension(0:), contiguous :: coefs_pb
!
!    coefs_pb(0) = coefs_bb(0)
!    coefs_pb(1) = - 3.0_PREC * coefs_bb(0) + 3.0_PREC * coefs_bb(1)
!    coefs_pb(2) = 3.0_PREC * coefs_bb(0) - 6.0_PREC * coefs_bb(1) &
!        + 3.0_PREC * coefs_bb(2)
!    coefs_pb(3) = -coefs_bb(0) + 3.0_PREC * coefs_bb(1) &
!        - 3.0_PREC * coefs_bb(2) + coefs_bb(3)
!
!end subroutine
!
!
!
!pure subroutine __APPEND(bernstein2power_deg2,__PREC) (coefs_bb, coefs_pb)
!
!    integer, parameter :: PREC = __PREC
!    real (PREC), intent(in), dimension(0:), contiguous :: coefs_bb
!    real (PREC), intent(out), dimension(0:), contiguous :: coefs_pb
!
!    coefs_pb(0) = coefs_bb(0)
!    coefs_pb(1) = 2.0_PREC * (coefs_bb(1) - coefs_bb(0))
!    coefs_pb(2) = coefs_bb(0) - 2.0_PREC * coefs_bb(1) + coefs_bb(2)
!
!end subroutine

