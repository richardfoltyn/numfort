module numfort_interpolate_fitpack08

    use iso_fortran_env
    ! interface specifications for F77 routines
    use fitpack_interfaces, only: splev_if, curfit_if
    use numfort_interpolate_result
    use numfort_common, only: workspace, ENUM_KIND
    implicit none
    private

    integer, parameter :: PREC = real64

    integer, parameter :: MIN_SPLINE_DEGREE = 1, MAX_SPLINE_DEGREE = 5, &
        DEFAULT_SPLINE_DEGREE = 3

    integer (ENUM_KIND), parameter :: SPLEV_EXTRAPOLATE = 0
    integer (ENUM_KIND), parameter :: SPLEV_ZERO = 1
    integer (ENUM_KIND), parameter :: SPLEV_BOUNDARY = 3

    public :: SPLEV_EXTRAPOLATE, SPLEV_ZERO, SPLEV_BOUNDARY

    interface splev
        module procedure splev_wrapper, splev_scalar
    end interface

    interface curfit
        module procedure curfit_wrapper
    end interface

    public :: curfit_get_nest, curfit, splev

contains

subroutine curfit_check_input (iopt, x, y, w, k, knots, coefs, stat, msg)

    integer, intent(in) :: iopt, k
    real (PREC), intent(in), dimension(:) :: x, y, w, knots, coefs
    integer, intent(out) :: stat
    character (len=*), intent(out) :: msg
    optional :: w

    integer :: nwrk, nest, m
    integer, parameter :: SUCCESS = 0

    stat = SUCCESS

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

    m = size(x)
    call curfit_get_wrk_size(m, k, nwrk, nest)

    if (size(knots) < nest .or. size(coefs) < nest) then
        msg = "knots or coefs array is too small"
        goto 100
    end if

    return

100 stat = INTERP_STATUS_INVALID_INPUT

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
    optional :: iopt, w, xb, xe, k, s, status, ssr, work
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

    call curfit_check_input (iopt, x, y, w, k, knots, coefs, istatus, msg)

    if (istatus == 0) then

        call curfit_get_wrk_size (m, k, nwrk, nest)

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
            knots, coefs, lssr, ptr_work%rwrk, nwrk, ptr_work%iwrk, lstatus)
    else
        print *, msg
    end if

    if (present(status)) status = lstatus
    if (present(ssr)) ssr = lssr

end subroutine

! ******************************************************************************
! CURFIT evaluation

subroutine splev_wrapper (knots, coefs, k, x, fx, ext, status)
    integer, intent(in) :: k
    real (PREC), intent(in), dimension(:), contiguous :: knots, coefs, x
    real (PREC), intent(out), dimension(:), contiguous :: fx
    integer, intent(in) :: ext
    integer, intent(out) :: status

    optional :: ext, status, k

    procedure (splev_if) :: splev

    integer :: m, n, lstatus, lext, lk

    lstatus = INTERP_STATUS_INVALID_INPUT

    m = size(x)
    n = size(knots)
    ! assume cubic splines by default
    lk = DEFAULT_SPLINE_DEGREE

    if (size(fx) /= m .or. size(coefs) /= n) then
        if (present(status)) status = lstatus
        return
    end if

    ! by default set function values outside of domain to boundary values
    lext = SPLEV_BOUNDARY

    if (present(k)) lk = k
    if (present(ext)) lext = ext

    call splev (knots, n, coefs, lk, x, fx, m, lext, lstatus)

    if (present(status)) status = lstatus

end subroutine

subroutine splev_scalar (knots, coefs, k, x, fx, ext, status)
    integer, intent(in) :: k
    real (PREC), intent(in), dimension(:), contiguous :: knots, coefs
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    integer, intent(in) :: ext
    integer, intent(out) :: status

    optional :: ext, status, k

    real (PREC) :: lx(1), lfx(1)

    lx(1) = x

    call splev_wrapper (knots, coefs, k, lx, lfx, ext, status)

    ! write back result
    fx = lfx(1)
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
