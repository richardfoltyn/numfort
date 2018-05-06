
pure subroutine __APPEND(pchip_fit_input_check,__PREC) (x, y, coef, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(in), dimension(:) :: coef
    type (status_t), intent(out) :: status

    integer :: n

    status = NF_STATUS_OK

    if (size(x) /= size(y)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    n = interp_pchip_get_ncoef (size(x))
    if (size(coef) < n) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

100 continue

end subroutine



pure subroutine __APPEND(interp_pchip_fit,__PREC) (x, y, coef, work, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(out), dimension(:) :: coef
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
    type (status_t), intent(out), optional :: status

    type (__APPEND(workspace,__PREC)), pointer :: ptr_work
    real (PREC), dimension(:), pointer, contiguous :: h, d, delta
    type (status_t) :: lstatus

    integer :: i, j, n, nrwrk
    real (PREC) :: h0, h1, delta0, delta1, w0, w1, d0, d1, ci, bi, sgn
    integer, parameter :: IDX_Y = 1, IDX_D = 2, IDX_C = 3, IDX_B = 4
        ! positions of individual coefficients within each interval block

    nullify (ptr_work, d, delta, h)

    lstatus = NF_STATUS_OK
    call pchip_fit_input_check (x, y, coef, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    n = size(x)
    ! Need arrays for h, delta, d
    nrwrk = n + 2 * (n-1)

    ! Allocate workspace arrays
    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal state in workspace object, in particular index offsets
    ! (this does not deallocate working arrays)
    call workspace_reset (ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk)

    call workspace_get_ptr (ptr_work, n-1, h)
    call workspace_get_ptr (ptr_work, n-1, delta)
    call workspace_get_ptr (ptr_work, n, d)

    ! Compute diff(x) and secant slopes for all intervals
    do i = 1, n-1
        h1 = x(i+1) - x(i)
        h(i) = h1
        delta(i) = (y(i+1) - y(i)) / h1
    end do

    ! Compute slopes at interior points
    do i = 2, n-1
        delta0 = delta(i-1)
        delta1 = delta(i)
        h0 = h(i-1)
        h1 = h(i)

        ! Evaluates to 1 of both have equal sign and are non-zero
        sgn = signum (delta0) * signum (delta1)
        w0 = 2.0_PREC * h1 + h0
        w1 = h1 + 2.0_PREC * h0
        d(i) = (w0 + w1) / (w0 / delta0) + w1/delta1
    end do

    ! Compute d(i) at endpoints
    d(1) = pchip_slope_end (h(1), h(2), delta(1), delta(2))
    d(n) = pchip_slope_end (h(n-2), h(n-1), delta(n-2), delta(n-1))

    ! Compute remaining coefficients C_i, B_i
    do i = 1, n-1
        delta0 = delta(i)
        h0 = h(i)
        d0 = d(i)
        d1 = d(i+1)

        ci = (3.0_PREC * delta0 - 2.0_PREC * d0 - d1) / h0
        bi = (d0 - 2.0_PREC * delta0 + d1) / h0**2.0_PREC

        ! Store everything in coefficient array
        j = (i-1) * (POLY_DEGREE+1)
        coef(j+IDX_Y) = y(i)
        coef(j+IDX_D) = d0
        coef(j+IDX_C) = ci
        coef(j+IDX_B) = bi
    end do

    ! store function value at endpoint for extrapolation
    coef(size(coef)) = y(n)

100 continue

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

    if (present(status)) status = lstatus

end subroutine


pure function __APPEND(pchip_slope_end,__PREC) (h0, h1, delta0, delta1) result(res)
    !*  PCHIP_SLOPE_END computes the slope at end points using the
    !   three-point formula
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: h0, h1
    real (PREC), intent(in) :: delta0, delta1
    real (PREC) :: res

    real (PREC) :: d

    d = ((2.0_PREC*h0 + h1)*delta0 - h0*delta1) / (h0 + h1)
    if (signum (d) /= signum(delta0)) then
        d = 0.0_PREC
    else if (signum(delta0) /= signum(delta1) .and. abs(d) > 3.0_PREC*abs(delta0)) then
        d = 3.0_PREC * delta0
    end if

    res = d
end function


pure subroutine __APPEND(pchip_eval_input_check,__PREC) (x, y, order, &
        ext, left, right, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: y
    integer, intent(in), optional :: order
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (status_t), intent(out) :: status

    integer, parameter :: ORDER_MAX = 2
    status = NF_STATUS_OK

    if (size(x) /= size(y)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(order)) then
        if (order < 0 .or. order > ORDER_MAX) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
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


pure subroutine __APPEND(interp_pchip_eval,__PREC) (xp, coef, x, y, order, ext, &
        left, right, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: coef
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: xp
        !*  x-values of data points. These must be the same points that were
        !   previously passed to INTERP_PCHIP_FIT.
    real (PREC), intent(out), dimension(:) :: y
    integer, intent(in), optional :: order
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left
    real (PREC), intent(in), optional :: right
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: lorder
    integer (NF_ENUM_KIND) :: lext
    real (PREC) :: lb, ub, yi, si, a
    integer :: i, j, jj, k, n

    real (PREC), parameter :: ALPHA(0:POLY_DEGREE,0:MAX_ORDER) = reshape( &
        [1, 1, 1, 1, 0, 1, 2, 3, 0, 0, 2, 6], &
        shape=[POLY_DEGREE+1, MAX_ORDER+1])

    lstatus = NF_STATUS_OK

    call pchip_eval_input_check (x, y, order, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lorder = 0
    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(order)) lorder = 0
    if (present(ext)) lext = ext

    lb = coef(1)
    ub = coef(size(coef))

    n = size(x)

    do i = 1, n
        if (x(i) < lb) then
            select case (lext)
            case (NF_INTERP_EVAL_BOUNDARY)
                y(i) = lb
                cycle
            case (NF_INTERP_EVAL_CONST)
                y(i) = left
                cycle
            case (NF_INTERP_EVAL_ERROR)
                lstatus = NF_STATUS_BOUNDS_ERROR
                goto 100
            end select
        else if (x(i) > 0) then
            select case (lext)
            case (NF_INTERP_EVAL_BOUNDARY)
                y(i) = ub
                cycle
            case (NF_INTERP_EVAL_CONST)
                y(i) = right
                cycle
            case (NF_INTERP_EVAL_ERROR)
                lstatus = NF_STATUS_BOUNDS_ERROR
                goto 100
            end select
        end if

        ! At this point we either have interior point or extrapolation
        ! was requested.
        ! Find interval j such that xp(j) <= x(i)
        j = interp_find (x(i), xp)
        jj = (j-1) * (POLY_DEGREE+1) + 1

        yi = 0.0
        si = x(i) - xp(j)
        do k = lorder, POLY_DEGREE
            a = ALPHA(k, lorder)
            yi = yi + a * coef(jj+k) * si ** (k-lorder)
        end do

    end do


100 continue

    if (present(status)) status = lstatus

end subroutine