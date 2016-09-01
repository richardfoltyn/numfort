module numfort_interpolate_fitpack08

    use iso_fortran_env
    use numfort_common, only: workspace, ENUM_KIND
    use numfort_interpolate_common
    use numfort_fitpack_interfaces
    use numfort_interpolate_result
    implicit none
    private

    integer, parameter :: PREC = real64

    integer, parameter :: MIN_SPLINE_DEGREE = 1, MAX_SPLINE_DEGREE = 5, &
        DEFAULT_SPLINE_DEGREE = 3

    ! used for internal input checking
    integer, parameter :: STATUS_INPUT_VALID = 0

    interface splev
        module procedure splev_wrapper, splev_scalar
    end interface

    interface curfit
        module procedure curfit_wrapper
    end interface

    interface splder
        module procedure splder_wrapper, splder_scalar
    end interface

    public :: curfit_get_nest, curfit, splev, splder

contains

! ******************************************************************************
! INPUT CHECKING used by various routines.

subroutine check_input (x, y, w, k, knots, coefs, stat, msg)

    integer, intent(in) :: k
    real (PREC), intent(in), dimension(:) :: x, y, w, knots, coefs
    integer, intent(out) :: stat
    character (len=*), intent(out) :: msg
    optional :: w

    stat = STATUS_INPUT_VALID

    if (size(x) /= size(y)) then
        msg = "Non-conformable arrays x, y"
        goto 100
    end if

    if (present(w)) then
        if (size(x) /= size(w)) then
            msg = "Non-conformable arrays x, w"
            goto 100
        end if
    end if

    ! check for supported spline degree
    if (k < MIN_SPLINE_DEGREE .or. k > MAX_SPLINE_DEGREE) then
        msg = "Unsupported spline degree k"
        goto 100
    end if

    if (size(knots) /= size(coefs)) then
        msg = "knots and coefs arrays must be of equal size"
        goto 100
    end if

    return

100 stat = INTERP_STATUS_INVALID_INPUT

end subroutine

subroutine check_input_ext (ext, stat, msg)
    integer :: ext, stat
    character (len=*) :: msg

    intent(in) :: ext
    intent(out) :: stat, msg

    integer :: i
    logical :: is_valid
    integer, parameter :: valid(4) = [INTERP_EVAL_ZERO, INTERP_EVAL_EXTRAPOLATE, &
        INTERP_EVAL_ERROR, INTERP_EVAL_BOUNDARY]

    is_valid = .false.

    do i = 1, size(valid)
        is_valid = is_valid .or. (ext == valid(i))
    end do

    if (.not. is_valid) then
        stat = INTERP_STATUS_INVALID_INPUT
        msg = "Invalid value extrapolation mode"
    else
        stat = STATUS_INPUT_VALID
    end if
end subroutine

! ******************************************************************************
! CURFIT fitting procedures

subroutine curfit_wrapper (iopt, x, y, w, xb, xe, k, s, work, n, knots, coefs, ssr, status)
    integer :: iopt, k, n, status
    real (PREC), dimension(:), contiguous :: x, y, w, knots, coefs
    real (PREC) :: s, xe, xb, ssr
    class (workspace), intent(in out), optional, target :: work

    intent(in) :: iopt, k, x, y, w, s, xe, xb
    intent(in out) :: knots, coefs, n, status, ssr
    optional :: iopt, w, xb, xe, k, s, status, ssr
    target :: w

    real (PREC), dimension(:), pointer, contiguous :: ptr_w
    integer :: m, liopt, lstatus, lk, nwrk, nest, istatus
    character (len=100) :: msg
    real (PREC) :: lxb, lxe, ls, lssr
    type (workspace), target :: lwork
    class (workspace), pointer :: ptr_work => null()

    procedure (curfit_if) :: curfit

    ! check for conformable arrays
    lstatus = INTERP_STATUS_INVALID_INPUT

    ! initialize default values
    lk = DEFAULT_SPLINE_DEGREE
    ! by default to LS on given set of knots
    liopt = -1
    m = size(x)
    lxb = x(1)
    lxe = x(m)
    ls = 0.0_PREC
    ! lower bound of recommended interval for s
    if (present(w)) ls = m - sqrt(2*real(m, PREC))

    ! override default values if any of the optional parameters were provided
    if (present(iopt)) liopt = iopt
    if (present(xb)) lxb = xb
    if (present(xe)) lxe = xe
    if (present(k)) lk = k
    if (present(s)) ls = s

    call check_input (x, y, w, lk, knots, coefs, istatus, msg)
    if (istatus /= STATUS_INPUT_VALID) goto 100

    call curfit_get_wrk_size (m, k, nwrk, nest)

    if (size(knots) < nest .or. size(coefs) < nest) then
        msg = "knots or coefs array is too small"
        goto 100
    end if

    if (present(work)) then
        ptr_work => work
    else
        ptr_work => lwork
    end if

    ! assert that working arrays are sufficiently large
    if (present(w)) then
        call ptr_work%assert_allocated (nrwrk=nwrk, niwrk=nest)
        ptr_w => w
    else
        ! allocate real working array such that we can append uniform
        ! weights at the end
        call ptr_work%assert_allocated (nrwrk=nwrk + m, niwrk=nest)
        ptr_w => ptr_work%rwrk(nwrk + 1:)
        ptr_w = 1.0_PREC
    end if

    ! call fitpack routine wrapper
    call curfit (liopt, m, x, y, ptr_w, lxb, lxe, lk, ls, nest, n, &
        knots, coefs, lssr, ptr_work%rwrk(1:nwrk), nwrk, ptr_work%iwrk, lstatus)

    if (present(status)) status = lstatus
    if (present(ssr)) ssr = lssr

    return

100 if(present(status)) status = INTERP_STATUS_INVALID_INPUT
    write (ERROR_UNIT, *) msg

end subroutine

! ******************************************************************************
! CURFIT evaluation

subroutine splev_wrapper (knots, coefs, k, x, y, ext, status)
    integer, intent(in) :: k
    real (PREC), intent(in), dimension(:), contiguous :: knots, coefs, x
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer, intent(in) :: ext
    integer, intent(out) :: status

    optional :: ext, status, k

    procedure (splev_if) :: splev

    integer :: m, n, lstatus, lext, lk, istatus
    character (100) :: msg

    lstatus = INTERP_STATUS_INVALID_INPUT

    m = size(x)
    n = size(knots)
    ! assume cubic splines by default
    lk = DEFAULT_SPLINE_DEGREE

    if (size(y) /= m) then
        msg = 'x and y array lengths differ'
        goto 100
    end if

    if (size(coefs) /= n) then
        msg = 'knots and coefs array lengths differ'
        goto 100
    end if

    ! by default set function values outside of domain to boundary values
    lext = INTERP_EVAL_EXTRAPOLATE

    if (present(k)) lk = k
    if (present(ext)) then
        call check_input_ext (ext, istatus, msg)
        if (istatus /= STATUS_INPUT_VALID) goto 100
        lext = ext
    end if

    call splev (knots, n, coefs, lk, x, y, m, lext, lstatus)

    if (present(status)) status = lstatus

    return

100 if (present(status)) status = INTERP_STATUS_INVALID_INPUT
    write (ERROR_UNIT, *) msg

end subroutine

subroutine splev_scalar (knots, coefs, k, x, y, ext, status)
    integer, intent(in) :: k
    real (PREC), intent(in), dimension(:), contiguous :: knots, coefs
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in) :: ext
    integer, intent(out) :: status

    optional :: ext, status, k

    real (PREC) :: lx(1), ly(1)

    lx(1) = x

    call splev_wrapper (knots, coefs, k, lx, ly, ext, status)

    ! write back result
    y = ly(1)
end subroutine

! ******************************************************************************
! Evaluation of derivatives

subroutine splder_wrapper (knots, coefs, k, order, x, y, ext, work, status)
    integer :: k, ext, status, order
    real (PREC), dimension(:), contiguous :: knots, coefs, x, y
    class (workspace), intent(in out), optional, target :: work

    intent (in) :: k, ext, knots, coefs, x, order
    intent (out) :: y, status

    optional :: ext, status, k, order

    type (workspace), target :: lwork
    type (workspace), pointer :: ptr_work

    integer :: lext, lorder, lstatus, lk, n, istatus, m
    character (100) :: msg

    procedure (splder_if) :: splder

    lk = DEFAULT_SPLINE_DEGREE
    lorder = 1
    lext = INTERP_EVAL_EXTRAPOLATE
    n = size(knots)
    m = size(x)

    if (present(k)) lk = k

    call check_input (x, y, k=lk, knots=knots, coefs=coefs, stat=istatus, msg=msg)
    if (istatus /= STATUS_INPUT_VALID) goto 100

    if (present(order)) then
        if (order < 0 .or. order > k) then
            msg = 'Invalid value: 0 <= order <= k required'
            goto 100
        end if
        lorder = order
    end if

    if (present(ext)) then
        call check_input_ext (ext, istatus, msg)
        if (istatus /= 0) goto 100
        lext = ext
    end if

    if (present(work)) then
        ptr_work => work
    else
        ptr_work => lwork
    end if
    call ptr_work%assert_allocated (nrwrk = n)

    call splder (knots, n, coefs, lk, lorder, x, y, m, lext, ptr_work%rwrk, lstatus)

    if (present(status)) status = lstatus
    return

100 if (present(status)) status = INTERP_STATUS_INVALID_INPUT
    write (ERROR_UNIT, *) msg
end subroutine

! SPLDER_SCALAR evaluates the derivate of a spline at a single scalar point
! x and stores the result in scalar y.
! Implemented as a wrapper for SPLDER_WRAPPER that creates local size-1 arrays.
subroutine splder_scalar (knots, coefs, k, order, x, y, ext, work, status)
    integer :: k, ext, status, order
    real (PREC), dimension(:), contiguous :: knots, coefs
    real (PREC) :: x, y
    class (workspace), intent(in out), optional, target :: work

    intent (in) :: k, ext, knots, coefs, x, order
    intent (out) :: y, status

    real (PREC), dimension(1) :: lx, ly

    lx(1) = x
    ly(1) = y

    call splder_wrapper (knots, coefs, k, order, lx, ly, ext, work, status)

    y = ly(1)
end subroutine

! ******************************************************************************
! Workspace query functions

function curfit_get_nest (m, k) result(res)
    integer, intent(in) :: k, m
    integer :: res

    if (k < MIN_SPLINE_DEGREE .or. k > MAX_SPLINE_DEGREE) then
        error stop "Unsupported spline degree"
    end if

    res = max(m + k + 1, 2*k + 3)
end function

subroutine curfit_get_wrk_size (m, k, nwrk, nest)
    integer, intent(in) :: m, k
    integer, intent(out) :: nwrk, nest

    nest = curfit_get_nest (m, k)
    nwrk = m * (k + 1) + nest * (7 + 3 * k)
end subroutine

end module
