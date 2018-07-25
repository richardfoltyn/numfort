
pure subroutine __APPEND(interp_linear_impl,__PREC) (x, xp, fp, fx, wgt, ext, left, right)
    !*  Implements linear interpolation for scalar input/return value
    !   and no optional arguments.
    !   Should be called from wrapper routines doing the optional argument handling,
    !   etc.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in) :: x
        !*  x-coordinate of the interpolated value
    real (PREC), intent(in), dimension(:) :: xp
        !*  x-coordinates of data points in increasing order
    real (PREC), intent(in), dimension(:) :: fp
        !*  y-coordinates of the data points, same length as xp.
    real (PREC), intent(out) :: fx
        !*  Interpolated value
    real (PREC), intent(out) :: wgt
        !*  Weight placed on the lower bound of the bracketing interval
    integer (NF_ENUM_KIND), intent(in) :: ext
        !*  Defines behavious in case function should be evaluated at point x
        !   outside of the range specified by [xp(1), xp(size(xp))].
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
        !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
        !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
        !           are used, respectively.
        !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and
        !           fp(size(xp)) for all x > xp(size(xp)).
        !       - NF_INTERP_EVAL_ZERO: Returns zero for all x outside of
        !           admissible range.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in) :: left
        !*  Value to return if x < xp(1)
    real (PREC), intent(in) :: right
        !*  Value to return if x > xp(size(xp))

    integer (INTSIZE) :: ilb, iub
    integer :: np

    np = size(xp)

    if (x < xp(1)) then
        wgt = 1.0_PREC
        select case (ext)
        case (NF_INTERP_EVAL_CONST)
            fx = left
            return
        case (NF_INTERP_EVAL_BOUNDARY)
            fx = fp(1)
            return
        case (NF_INTERP_EVAL_ZERO)
            fx = 0.0_PREC
            return
        end select
    else if (x > xp(np)) then
        wgt = 0.0_PREC
        select case (ext)
        case (NF_INTERP_EVAL_CONST)
            fx = right
            return
        case (NF_INTERP_EVAL_BOUNDARY)
            fx = fp(np)
            return
        case (NF_INTERP_EVAL_ZERO)
            fx = 0.0_PREC
            return
        end select
    end if

    ! default: At this point there is either a bracketing interval, or we extrapolate
    ! values outside of domain
    ilb = interp_find (x, xp)
    iub = ilb + 1
    ! Weight on lower bound
    wgt = (xp(iub) - x) / (xp(iub) - xp(ilb))
    fx = wgt * fp(ilb) + (1.0_PREC-wgt) * fp(iub)

end subroutine



pure subroutine __APPEND(interp_linear_scalar,__PREC) (x, xp, fp, fx, ext, &
        left, right, wgt)
    !*  INTERP_LINEAR_SCALAR performs linear interpolation for scalar arguments
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in) :: x
        !*  The x-coordinate of the interpolated value.
    real (PREC), intent(in), dimension(:) :: xp
        !*  The x-coordinates of data points in increasing order.
    real (PREC), intent(in), dimension(:) :: fp
        !*  The y-coordinates of data points, same length as xp.
    real (PREC), intent(out) :: fx
        !*  Interpolated value
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point x
        !   outside of the range specified by [xp(1), xp(size(xp))].
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
        !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
        !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
        !           are used, respectively.
        !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and
        !           fp(size(xp)) for all x > xp(size(xp)).
        !       - NF_INTERP_EVAL_ZERO: Returns zero for all x outside of
        !           admissible range.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: left
        !*  Value to return if x < xp(1) and extrapolation is disabled (`ext=.false.`).
        !   Default is fp(1).
    real (PREC), intent(in), optional :: right
        !*  Value to return if x > xp(-1) and extrapolation disabled (`ext=.false.`).
        !   Default is fp(-1).
    real (PREC), intent(out), optional :: wgt

    integer (NF_ENUM_KIND) :: lext
    real (PREC) :: lright, lleft, lwgt

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    if (present(left)) then
        lleft = left
    else
        lleft = fp(1)
    end if

    if (present(right)) then
        lright = right
    else
        lright = fp(size(fp))
    end if

    call interp_linear_impl (x, xp, fp, fx, lwgt, lext, lleft, lright)

    if (present(wgt)) wgt = lwgt

end subroutine



pure subroutine __APPEND(interp_linear_1d,__PREC) (x, xp, fp, fx, ext, left, right, wgt)
    !*  INTERP_LINEAR_1D performs linear interpolation for an array of
    !   points.

    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: x
        !*  The x-coordinates of the interpolated values.
    real (PREC), intent(in), dimension(:) :: xp
        !*  The x-coordinates of data points in increasing order.
    real (PREC), intent(in), dimension(:) :: fp
        !*  The y-coordinates of data points, same length as xp.
    real (PREC), intent(out), dimension(:) :: fx
        !*  Array of interpolated values, same shape as x.
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point x
        !   outside of the range specified by [xp(1), xp(size(xp))].
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
        !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
        !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and
        !           fp(size(xp)) for all x > xp(size(xp)).
        !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
        !           are used, respectively.
        !       - NF_INTERP_EVAL_ZERO: Returns zero for all x outside of
        !           admissible range.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: left
        !*  Value to return if x < xp(1) and extrapolation is disabled (`ext=.false.`).
        !   Default is fp(1)
    real (PREC), intent(in), optional :: right
        !*  Value to return if x > xp(-1) and extrapolation disabled (`ext=.false.`).
        !   Default is fp(-1)
    real (PREC), intent(out), dimension(:), optional :: wgt

    integer (INTSIZE) :: n, i
    integer (NF_ENUM_KIND) :: lext
    real (PREC) :: lright, lleft, wgt_dummy

    n = size (x)
    lext = NF_INTERP_EVAL_EXTRAPOLATE

    if (present(ext)) lext = ext
    if (present(left)) then
        lleft = left
    else
        lleft = fp(1)
    end if

    if (present(right)) then
        lright = right
    else
        lright = fp(size(fp))
    end if

    if (present(wgt)) then
        do i = 1, n
            ! call scalar implementation
            call interp_linear_impl (x(i), xp, fp, fx(i), wgt(i), lext, &
                lleft, lright)
        end do
    else
        do i = 1, n
            ! call scalar implementation
            call interp_linear_impl (x(i), xp, fp, fx(i), wgt_dummy, lext, &
                lleft, lright)
        end do
    end if

end subroutine


pure subroutine __APPEND(interp_bilinear_impl,__PREC) (x1, x2, xp1, xp2, fp, &
        ext, fill_value, fx, wgt)
    !*  INTERP_BILINEAR_IMPL performs bilineare interpolation over a
    !   rectangular grid.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in) :: x1
    real (PREC), intent(in) :: x2
    real (PREC), intent(in), dimension(:) :: xp1
        !*  Grid points in dimension 1
    real (PREC), intent(in), dimension(:) :: xp2
        !*  Grid points in dimension 2
    real (PREC), intent(in), dimension(:,:) :: fp
        !*  Function evaluated at given grid points, ie FP(i,j) = f(xp1(i),xp2(j))
    integer (NF_ENUM_KIND), intent(in) :: ext
        !*  Defines behavious in case function should be evaluated at point (x1,x2)
        !   outside of the range specified by [xpi(1), xpi(size(xpi))] for
        !   i = 1,2.
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'fill_value'
        !           for all xi < xpi(1), i=1,2.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in) :: fill_value
        !*  Fill value to be used if point (X1,X2) is outside of the domain
        !   defined by XP1, XP2.
    real (PREC), intent(out) :: fx
        !*  Interpolated function value f(x1,x2)
    real (PREC), intent(out), dimension(:) :: wgt
        !*  If present, returns the weights for the lower bound of the bracketing
        !   interval for each dimension.

    integer (INTSIZE) :: ilb1, ilb2
    real (PREC) :: w1, w2, fx1_lb, fx1_ub
    integer :: np1, np2

    np1 = size(xp1)
    np2 = size(xp2)

    select case (ext)
    case (NF_INTERP_EVAL_CONST)
        if (x1 < xp1(1) .or. x1 > xp1(np1) .or. x2 < xp2(1) .or. x2 > xp2(np2)) then
            wgt = 0.0_PREC
            goto 100
        end if
    end select

    ! default: At this point there is either a bracketing interval, or we extrapolate
    ! values outside of domain
    ilb1 = interp_find (x1, xp1)
    ilb2 = interp_find (x2, xp2)

    w1 = (x1 - xp1(ilb1)) / (xp1(ilb1+1) - xp1(ilb1))
    w2 = (x2 - xp2(ilb2)) / (xp2(ilb2+1) - xp2(ilb2))

    ! Interpolate in dimension 1 for upper and lower bracketing value of
    ! dimension 2.
    fx1_lb = (1.0_PREC-w1) * fp(ilb1,ilb2) + w1 * fp(ilb1+1,ilb2)
    fx1_ub = (1.0_PREC-w1) * fp(ilb1,ilb2+1) + w1 * fp(ilb1+1,ilb2+1)

    ! Interpolate in dimension two
    fx = (1.0_PREC-w2) * fx1_lb + w2 * fx1_ub

    ! Return weights on lower bounds
    wgt(1) = 1.0_PREC - w1
    wgt(2) = 1.0_PREC - w2

    return

100 continue
    ! Return constant fill_value if X in any dimension is not "interior".

    fx = fill_value

end subroutine



pure subroutine __APPEND(interp_bilinear_check_input,__PREC) (xp1, xp2, fp, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: xp1, xp2
    real (PREC), intent(in), dimension(:,:) :: fp
    type (status_t), intent(out) :: status

    status = NF_STATUS_OK

    if (size(xp1) /= size(fp, 1) .or. size(xp2) /= size(fp,2)) then
        status = NF_STATUS_INVALID_ARG
        return
    end if

    if (size(xp1) < 2 .or. size(xp2) < 2) then
        status = NF_STATUS_INVALID_ARG
        return
    end if

end subroutine



pure subroutine __APPEND(interp_bilinear_1d,__PREC) (x1, x2, xp1, xp2, fp, fx, &
        ext, fill_value, wgt, status)
    !*  INTERP_BILINEAR_IMPL performs bilineare interpolation over a
    !   rectangular grid.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:) :: x1
        !*  Array of dimension-1 values where function should be interpolated
    real (PREC), intent(in), dimension(:) :: x2
        !*  Array of dimension-2 values where function should be interpolated
    real (PREC), intent(in), dimension(:) :: xp1
        !*  Grid points in dimension 1
    real (PREC), intent(in), dimension(:) :: xp2
        !*  Grid points in dimension 2
    real (PREC), intent(in), dimension(:,:) :: fp
        !*  Function evaluated at given grid points, ie FP(i,j) = f(xp1(i),xp2(j))
    real (PREC), intent(out), dimension(:) :: fx
        !*  Array of interpolated values
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point (x1,x2)
        !   outside of the range specified by [xpi(1), xpi(size(xpi))] for
        !   i = 1,2.
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'fill_value'
        !           for all xi < xpi(1), i=1,2.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: fill_value
        !*  Fill value to be used if point (X1,X2) is outside of the domain
        !   defined by XP1, XP2.
    real (PREC), intent(out), dimension(:,:), optional :: wgt
    type (status_t), intent(out), optional :: status
        !*  Optional exit status.

    type (status_t) :: lstatus
    real (PREC) :: lfill_value
    integer (NF_ENUM_KIND) :: lext
    integer (INTSIZE) :: i
    real (PREC), dimension(2) :: wgt_dummy

    lstatus = NF_STATUS_OK

    if (size(x1) /= size(x2) .or. size(x1) /= size(fx)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Check remaining inputs using common function
    call interp_bilinear_check_input (xp1, xp2, fp, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(wgt)) then
        if (size(wgt, 1) < 2) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    lfill_value = 0.0_PREC
    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext
    if (present(fill_value)) then
        lfill_value = fill_value
        ! Assume this implies setting non-interior points to constant
        lext = NF_INTERP_EVAL_CONST
    end if

    if (present(wgt)) then
        do i = 1, size(x1)
            call interp_bilinear_impl (x1(i), x2(i), xp1, xp2, fp, lext, &
                lfill_value, fx(i), wgt(:,i))
        end do
    else
        do i = 1, size(x1)
            call interp_bilinear_impl (x1(i), x2(i), xp1, xp2, fp, lext, &
                lfill_value, fx(i), wgt_dummy)
        end do
    end if

100 continue
    if (present(status)) status = lstatus
end subroutine



pure subroutine __APPEND(interp_bilinear_scalar,__PREC) (x1, x2, xp1, xp2, fp, fx, &
        ext, fill_value, wgt, status)
    !*  INTERP_BILINEAR_IMPL performs bilineare interpolation over a
    !   rectangular grid.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in) :: x1
        !*  Dimension-1 value where function should be interpolated
    real (PREC), intent(in) :: x2
        !*  Dimension-2 value where function should be interpolated
    real (PREC), intent(in), dimension(:) :: xp1
        !*  Grid points in dimension 1
    real (PREC), intent(in), dimension(:) :: xp2
        !*  Grid points in dimension 2
    real (PREC), intent(in), dimension(:,:) :: fp
        !*  Function evaluated at given grid points, ie FP(i,j) = f(xp1(i),xp2(j))
    real (PREC), intent(out) :: fx
        !*  Array of interpolated values
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point (x1,x2)
        !   outside of the range specified by [xpi(1), xpi(size(xpi))] for
        !   i = 1,2.
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'fill_value'
        !           for all xi < xpi(1), i=1,2.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: fill_value
        !*  Fill value to be used if point (X1,X2) is outside of the domain
        !   defined by XP1, XP2.
    real (PREC), intent(out), dimension(:), optional :: wgt
    type (status_t), intent(out), optional :: status
        !*  Optional exit status.

    type (status_t) :: lstatus
    real (PREC) :: lfill_value
    integer (NF_ENUM_KIND) :: lext
    real (PREC), dimension(2) :: wgt_dummy

    lstatus = NF_STATUS_OK

    ! Check remaining inputs using common function
    call interp_bilinear_check_input (xp1, xp2, fp, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(wgt)) then
    if (size(wgt) < 2) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    end if
    lfill_value = 0.0_PREC
    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext
    if (present(fill_value)) then
        lfill_value = fill_value
        ! Assume this implies setting non-interior points to constant
        lext = NF_INTERP_EVAL_CONST
    end if

    if (present(wgt)) then
        call interp_bilinear_impl (x1, x2, xp1, xp2, fp, lext, lfill_value, &
            fx, wgt)
    else
        call interp_bilinear_impl (x1, x2, xp1, xp2, fp, lext, lfill_value, &
            fx, wgt_dummy)
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine


