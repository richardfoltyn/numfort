

pure subroutine __APPEND(ppoly2d_val_check_input,__PREC) (self, x1, x2, y, &
        ext, ext_val, status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: x1
    real (PREC), intent(in), dimension(:) :: x2
    real (PREC), intent(out), dimension(:) :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: ext_val
    type (status_t), intent(out) :: status

    integer (NF_ENUM_KIND), parameter :: EXT_VALID(*) = &
        [   NF_INTERP_EVAL_EXTRAPOLATE, &
            NF_INTERP_EVAL_ERROR, &
            NF_INTERP_EVAL_CONST ]

    status = NF_STATUS_oK

    if (size(x1) /= size(x2) .or. size(x1) /= size(y)) goto 100

    if (present(ext)) then
        if (all(ext /= EXT_VALID)) goto 100
        if (.not. present(ext_val)) goto 100
    end if

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



pure subroutine __APPEND(ppoly2d_val,__PREC) (self, knots, coefs, x1, x2, y, &
        ext, ext_val, status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(in) :: self
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    real (PREC), intent(in), dimension(:) :: x1
    real (PREC), intent(in), dimension(:) :: x2
    real (PREC), intent(out), dimension(:) :: y
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: ext_val
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer (NF_ENUM_KIND) :: lext
    real (PREC) :: lext_val

    call ppoly2d_check_input (self, self%nknots, self%degree, knots, coefs, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call ppoly2d_val_check_input (self, x1, x2, y, ext, ext_val, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    lext_val = 0.0
    if (present(ext)) lext = ext
    if (present(ext_val)) lext_val = ext_val

    select case (self%degree)
    case (1)
        call ppoly2d_val_bilinear (self, knots, coefs, x1, x2, y, lext, lext_val, lstatus)
    case default
        lstatus = NF_STATUS_NOT_IMPLEMENTED
    end select

100 continue

    if (present(status)) status = lstatus
end subroutine


pure subroutine __APPEND(ppoly2d_val_bilinear,__PREC) (self, knots, coefs, x1, x2, y, &
        ext, ext_val, status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(in) :: self
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
    real (PREC), intent(in), dimension(:) :: x1
    real (PREC), intent(in), dimension(:) :: x2
    real (PREC), intent(out), dimension(:) :: y
    integer (NF_ENUM_KIND), intent(in) :: ext
    real (PREC), intent(in) :: ext_val
    type (status_t), intent(out) :: status

    integer :: nx, n1, n2, k
    integer :: i1, i2, i, jj
    real (PREC) :: x1i, x2i, x1lb, x2lb, x1ub, x2ub, z1, z2

    nx = size(x1)
    n1 = self%nknots(1)
    n2 = self%nknots(2)
    k = self%degree

    x1lb = knots(1)
    x1ub = knots(n1)
    x2lb = knots(n1+1)
    x2ub = knots(n1+n2)

    do i = 1, nx
        x1i = x1(i)
        if (x1i < x1lb) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                i1 = 1
                x1i = x1lb
            case (NF_INTERP_EVAL_CONST)
                y(i) = ext_val
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval that will be used to interpolate/extrapolate
                ! in dimension 1
                i1 = interp_find (x1i, knots(1:n1))
            end select
        else if (x1i > x1ub) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                i1 = n1 - 1
                x1i = x1ub
            case (NF_INTERP_EVAL_CONST)
                y(i) = ext_val
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval that will be used to interpolate/extrapolate
                ! in dimension 1
                i1 = interp_find (x1i, knots(1:n1))
            end select
        else
            i1 = interp_find (x1i, knots(1:n1))
        end if

        x2i = x2(i)
        if (x2i < x2lb) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                i2 = 1
                x2i = x2lb
            case (NF_INTERP_EVAL_CONST)
                y(i) = ext_val
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval that will be used to interpolate/extrapolate
                ! in dimension 1
                i2 = interp_find (x2i, knots(n1+1:n1+n2))
            end select
        else if (x2i > x2ub) then
            select case (ext)
            case (NF_INTERP_EVAL_BOUNDARY)
                i2 = n1 - 1
                x2i = x2ub
            case (NF_INTERP_EVAL_CONST)
                y(i) = ext_val
                cycle
            case (NF_INTERP_EVAL_ERROR)
                status = NF_STATUS_BOUNDS_ERROR
                goto 100
            case default
                ! Find interval that will be used to interpolate/extrapolate
                ! in dimension 1
                i2 = interp_find (x2i, knots(n1+1:n1+n2))
            end select
        else
            i2 = interp_find (x2i, knots(n1+1:n1+n2))
        end if

        ! Offset of coefficient block for bracket (i1,i2)
        jj = ((i2-1) * n1 + (i1-1)) * (k+1)**2

        z1 = (x1i - knots(i1)) / (knots(i1+1) - knots(i1))
        z2 = (x2i - knots(n1+i2)) / (knots(n1+i2+1) - knots(n1+i2))

        y(i) = coefs(jj) + coefs(jj+1)*z1 + coefs(jj+2)*z2 + coefs(jj+3)*z1*z2

    end do

100 continue

end subroutine

