
subroutine __APPEND(interp_linear_impl,__PREC) (x, xp, fp, fx, ext, left, right)
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
    real (PREC), intent(out) :: fx
        !*  Interpolated value

    integer (INTSIZE) :: ilb, iub
    real (PREC) :: slope
    integer :: np

    np = size(xp)

    if (x < xp(1)) then
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
    slope = (x - xp(ilb)) / (xp(iub) - xp(ilb))
    fx = (1.0_PREC-slope) * fp(ilb) + slope * fp(iub)

end subroutine



subroutine __APPEND(interp_linear_scalar,__PREC) (x, xp, fp, fx, ext, left, right)
    !*  INTERP_LINEAR_SCALAR performs linear interpolation for scalar arguments
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in) :: x
        !*  The x-coordinate of the interpolated value.
    real (PREC), intent(in), dimension(:) :: xp
        !*  The x-coordinates of data points in increasing order.
    real (PREC), intent(in), dimension(:) :: fp
        !*  The y-coordinates of data points, same length as xp.
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
    real (PREC), intent(out) :: fx
        !*  Interpolated value

    integer (NF_ENUM_KIND) :: lext
    real (PREC) :: lright, lleft

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

    call interp_linear_impl (x, xp, fp, fx, lext, lleft, lright)

end subroutine



subroutine __APPEND(interp_linear_1d,__PREC) (x, xp, fp, fx, ext, left, right)
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
    real (PREC), intent(out), dimension(:) :: fx
        !*  Array of interpolated values, same shape as x.

    integer (INTSIZE) :: n, i
    integer (NF_ENUM_KIND) :: lext
    real (PREC) :: lright, lleft

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

    do i = 1, n
        ! call scalar implementation
        call interp_linear_impl (x(i), xp, fp, fx(i), lext, lleft, lright)
    end do

end subroutine


