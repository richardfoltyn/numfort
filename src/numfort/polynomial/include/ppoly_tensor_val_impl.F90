

pure subroutine __APPEND(ppoly2d_val_check_input,__PREC) (self, x1, x2, y, &
        ext, ext_value, m, dim, status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: x1
    real (PREC), intent(in), dimension(:) :: x2
    real (PREC), intent(in), dimension(:) :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: ext_value
    integer, intent(in), optional :: m
    integer, intent(in), optional :: dim
    type (status_t), intent(out) :: status

    integer (NF_ENUM_KIND), parameter :: EXT_VALID(*) = &
        [   NF_INTERP_EVAL_EXTRAPOLATE, &
            NF_INTERP_EVAL_ERROR, &
            NF_INTERP_EVAL_CONST ]

    status = NF_STATUS_oK

    if (size(x1) /= size(x2) .or. size(x1) /= size(y)) goto 100

    if (present(m)) then
        if (.not. present(dim)) goto 100
        ! Note: at this point we know that DIM is present
        if (dim /= 1 .and. dim /= 2) goto 100
        if (m < 0) goto 100
    end if

    if (present(ext)) then
        if (all(ext /= EXT_VALID)) goto 100
        if ((ext == NF_INTERP_EVAL_CONST) .and. .not. present(ext_value)) goto 100
    end if

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



pure subroutine __APPEND(ppoly2d_val,__PREC) (self, knots, coefs, x1, x2, y, &
        cache, ext, ext_value, m, dim, status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(in) :: self
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    real (PREC), intent(in), dimension(:) :: x1
    real (PREC), intent(in), dimension(:) :: x2
    real (PREC), intent(out), dimension(:) :: y
    type (search_cache), intent(inout), dimension(:), optional :: cache
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: ext_value
    integer, intent(in), optional :: m
    integer, intent(in), optional :: dim
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer (NF_ENUM_KIND) :: lext
    real (PREC) :: lext_value
    type (search_cache), dimension(2) :: lcache
    integer :: lm, ldim

    call ppoly2d_check_input (self, self%nknots, self%degree, knots, coefs, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call ppoly2d_val_check_input (self, x1, x2, y, ext, ext_value, &
        m, dim, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    lext_value = 0.0
    lm = 0
    ldim = 0
    if (present(ext)) lext = ext
    if (present(ext_value)) lext_value = ext_value
    if (present(cache)) lcache = cache
    if (present(m)) lm = m
    if (present(dim)) ldim = dim

    if (lm > 1) then
        ! TODO: implement support for higher-order polynomials
        lstatus = NF_STATUS_NOT_IMPLEMENTED
        goto 100
    end if

    select case (self%degree)
    case (1)
        call ppoly2d_val_bilinear (self, knots, coefs, x1, x2, y, lcache, &
            lext, lext_value, lm, ldim, lstatus)
    case default
        lstatus = NF_STATUS_NOT_IMPLEMENTED
    end select

100 continue

    if (present(status)) status = lstatus
    if (present(cache)) cache = lcache
end subroutine



pure subroutine __APPEND(ppoly2d_val_bilinear,__PREC) (self, knots, coefs, &
        x1, x2, y, cache, ext, ext_value, m, dim, status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(in) :: self
    real (PREC), intent(in), dimension(:), contiguous :: knots
        !*  Knots defining the breakpoints of the piecewise bivariate polynomial
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
        !*  Piecewise bivariate polynomial coefficients
    real (PREC), intent(in), dimension(:) :: x1
        !*  x1-coordinates of points where function should be interpolated
    real (PREC), intent(in), dimension(:) :: x2
        !*  x2-coordinates of points where function should be interpolated
    real (PREC), intent(out), dimension(:) :: y
        !*  Interpolated function value (or its derivative)
    type (search_cache), intent(inout), dimension(:) :: cache
        !*  Optional search cache to accelerate repeated searches
    integer (NF_ENUM_KIND), intent(in) :: ext
        !*  Flag controlling extrapolation. Valid values are
        !       NF_INTERP_EVAL_ERROR:
        !           Return NF_STATUS_BOUNDS_ERROR whenever a non-interior
        !           value (x1,x2) is encountered.
        !       NF_INTERP_EVAL_CONST:
        !           Set Y(i) = EXT_VALUE for each i with non-interior
        !           coordinates (x1(i),x2(i))
        !       NF_INTERP_EVAL_EXTRAPOLATE:
        !           Extrapolate value using the bilinear polynomial in the
        !           most adjacent "bounding" rectangle.
    real (PREC), intent(in) :: ext_value
        !*  Value assigned to non-interior points when EXT=NF_INTER_EVAL_CONST.
    integer, intent(in) :: m
        !*  (partial) derivative order. M=0 computes the function value.
    integer, intent(in) :: dim
        !*  Dimension in which partial derivative is taken.
    type (status_t), intent(out) :: status
        !*  Exit status

    integer :: nx, n1, n2
    integer :: i1, i2, i, jj, j
    real (PREC) :: x1i, x2i, x1lb, x2lb, x1ub, x2ub, z1, z2, dx1, dx2, yi
    integer, parameter :: NCOEFS_BLOCK = 4
    real (PREC), dimension(0:NCOEFS_BLOCK-1) :: zz

    status = NF_STATUS_OK

    nx = size(x1)
    n1 = self%nknots(1)
    n2 = self%nknots(2)

    x1lb = knots(1)
    x1ub = knots(n1)
    x2lb = knots(n1+1)
    x2ub = knots(n1+n2)

    zz(0) = 1.0_PREC

    do i = 1, nx
        x1i = x1(i)
        if (x1i < x1lb) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                i1 = 1
                x1i = x1lb
            case (NF_INTERP_EVAL_CONST)
                y(i) = ext_value
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval that will be used to interpolate/extrapolate
                ! in dimension 1
                call bsearch_cached (x1i, knots(1:n1), i1, cache(1))
            end select
        else if (x1i > x1ub) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                i1 = n1 - 1
                x1i = x1ub
            case (NF_INTERP_EVAL_CONST)
                y(i) = ext_value
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval that will be used to interpolate/extrapolate
                ! in dimension 1
                call bsearch_cached (x1i, knots(1:n1), i1, cache(1))
            end select
        else
            call bsearch_cached (x1i, knots(1:n1), i1, cache(1))
        end if

        x2i = x2(i)
        if (x2i < x2lb) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                i2 = 1
                x2i = x2lb
            case (NF_INTERP_EVAL_CONST)
                y(i) = ext_value
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval that will be used to interpolate/extrapolate
                ! in dimension 1
                call bsearch_cached (x2i, knots(n1+1:n1+n2), i2, cache(2))
            end select
        else if (x2i > x2ub) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                i2 = n1 - 1
                x2i = x2ub
            case (NF_INTERP_EVAL_CONST)
                y(i) = ext_value
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval that will be used to interpolate/extrapolate:
                ! in dimension 1
                call bsearch_cached (x2i, knots(n1+1:n1+n2), i2, cache(2))
            end select
        else
            call bsearch_cached (x2i, knots(n1+1:n1+n2), i2, cache(2))
        end if

        ! Offset of coefficient block for bracket (i1,i2)
        jj = ((i2-1) * (n1-1) + (i1-1)) * NCOEFS_BLOCK

        dx1 = knots(i1+1) - knots(i1)
        dx2 = knots(n1+i2+1) - knots(n1+i2)
        zz(1) = (x1i - knots(i1)) / dx1
        zz(2) = (x2i - knots(n1+i2)) / dx2

        if (m == 0) then
            zz(3) = zz(1)*zz(2)
            yi = coefs(jj)
            do j = 1, NCOEFS_BLOCK-1
                yi = yi + coefs(jj+j)*zz(j)
            end do
            y(i) = yi
        else if (m == 1) then
            if (dim == 1) then
                ! Apply chain rule dz1/dx1
                y(i) = (coefs(jj+1) + coefs(jj+3)*zz(2)) / dx1
            else if (dim == 2) then
                ! Apply chain rule dz2/dx2
                y(i) = (coefs(jj+2) + coefs(jj+3)*zz(1)) / dx2
            end if
        end if

    end do

100 continue

end subroutine



pure subroutine __APPEND(ppoly2d_val_scalar,__PREC) (self, knots, &
        coefs, x1, x2, y, cache, ext, ext_value, m, dim, status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(in) :: self
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
    real (PREC), intent(in) :: x1
    real (PREC), intent(in) :: x2
    real (PREC), intent(out) :: y
    type (search_cache), intent(inout), dimension(:), optional :: cache
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: ext_value
    integer, intent(in), optional :: m
    integer, intent(in), optional :: dim
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(1) :: x1a, x2a, ya

    x1a(1) = x1
    x2a(1) = x2
    ya(1) = 0.0

    call ppolyval (self, knots, coefs, x1a, x2a, ya, cache, ext, ext_value, &
        m, dim, status)
    y = ya(1)

end subroutine

