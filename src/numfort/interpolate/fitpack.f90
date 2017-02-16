module numfort_interpolate_fitpack

    use iso_fortran_env
    use numfort_common, only: workspace, ENUM_KIND
    use numfort_interpolate_common
    use numfort_interpolate_result
    use fitpack_real64, only: fitpack_curfit_real64 => curfit, &
        fitpack_concon_real64 => concon, fitpack_splev_real64 => splev, &
        fitpack_splder_real64 => splder

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

    interface concon
        module procedure concon_wrapper
    end interface

    interface splder
        module procedure splder_wrapper, splder_scalar
    end interface

    public :: curfit_get_nest, curfit, splev, splder
    public :: concon_get_nest, concon

contains

! ******************************************************************************
! INPUT CHECKING used by various routines.

!>  Verify that arrays ore of equal length, conditional on being present.
!   If arr2 is not present, array length is considered to be equal.
pure function check_length (arr1, arr2) result(status)
    real (PREC), intent(in), dimension(:) :: arr1
    real (PREC), intent(in), dimension(:), optional :: arr2
    integer (ENUM_KIND) :: status

    status = INTERP_STATUS_SUCCESS
    if (present(arr2)) then
        if (size(arr1) /= size(arr2)) then
            status = INTERP_STATUS_INVALID_INPUT
        end if
    end if

end function

pure subroutine check_input (x, y, w, k, knots, coefs, v, sx, bind, status, msg)

    integer, intent(in), optional :: k
    real (PREC), intent(in), dimension(:) :: x, y, knots, coefs
    real (PREC), intent(in), dimension(:), optional :: w
        !!  Vector of weights, defaults to 1.0
    real (PREC), intent(in), dimension(:), optional :: v
        !!  Concavity restrictions, used by concon
    real (PREC), intent(in), dimension(:), optional :: sx
        !!  Used by concon
    logical, intent(in), dimension(:), optional :: bind
        !!  Array containing flag whether constraint is binding: used for concon

    integer (ENUM_KIND), intent(out) :: status
    character (len=*), intent(out), optional :: msg

    status = INTERP_STATUS_SUCCESS

    if (size(x) /= size(y)) then
        if (present(msg)) msg = "Non-conformable arrays x, y"
        goto 100
    end if

    if (check_length (x, w) /= INTERP_STATUS_SUCCESS) then
        if (present(msg)) msg = "Non-conformable arrays x, w"
        goto 100
    end if

    ! check for supported spline degree
    if (present(k)) then
        if (k < MIN_SPLINE_DEGREE .or. k > MAX_SPLINE_DEGREE) then
            if (present(msg)) msg = "Unsupported spline degree k"
            goto 100
        end if
    end if

    if (size(knots) /= size(coefs)) then
        if (present(msg)) msg = "knots and coefs arrays must be of equal size"
        goto 100
    end if

    if (check_length (x, sx) /= INTERP_STATUS_SUCCESS) then
        if (present(msg)) msg = "x and sx must be of equal size"
        goto 100
    end if

    ! bind is logical array, don't bother to create check_length() that
    ! accepts logical arguments.
    if (present(bind)) then
        if (size(bind) /= size(knots)) then
            if (present(msg)) msg = "knots and bind arrays must be of equal size"
            goto 100
        end if
    end if

    if (check_length (x, v) /= INTERP_STATUS_SUCCESS) then
        if (present(msg)) msg = "x and v must be of equal size"
        goto 100
    end if

    return

100 status = INTERP_STATUS_INVALID_INPUT

end subroutine

pure subroutine check_input_ext (ext, status, msg)
    integer, intent(in) :: ext
    integer (ENUM_KIND), intent(out) :: status
    character (len=*), intent(out), optional :: msg

    integer :: i
    logical :: is_valid
    integer, parameter :: valid(4) = [INTERP_EVAL_ZERO, INTERP_EVAL_EXTRAPOLATE, &
        INTERP_EVAL_ERROR, INTERP_EVAL_BOUNDARY]

    is_valid = .false.

    do i = 1, size(valid)
        is_valid = is_valid .or. (ext == valid(i))
    end do

    if (.not. is_valid) then
        status = INTERP_STATUS_INVALID_INPUT
        if (present(msg)) msg = "Invalid value extrapolation mode"
    else
        status = STATUS_INPUT_VALID
    end if
end subroutine

! ******************************************************************************
! CURFIT fitting procedures

pure subroutine curfit_wrapper (iopt, x, y, w, xb, xe, k, s, work, n, knots, &
        coefs, ssr, status, msg)
    integer, intent(in), optional :: iopt
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: y
    real (PREC), intent(in), dimension(:), contiguous, optional :: w
    real (PREC), intent(in), optional :: xe
    real (PREC), intent(in), optional :: xb
    integer, intent(in), optional :: k
    real (PREC), intent(in), optional :: s
    class (workspace), intent(in out), optional, target :: work
    integer, intent(out) :: n
    real (PREC), intent(in out), dimension(:), contiguous :: knots
    real (PREC), intent(out), dimension(:), contiguous :: coefs
    real (PREC), intent(out), optional :: ssr
    integer (ENUM_KIND), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: m, liopt, lstatus, lk, nwrk, nest, istatus
    real (PREC) :: lxb, lxe, ls, lssr
    type (workspace), target :: lwork
    class (workspace), pointer :: ptr_work

    ! check for conformable arrays
    lstatus = INTERP_STATUS_INVALID_INPUT
    nullify(ptr_work)

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

    call check_input (x, y, w, lk, knots, coefs, status=istatus, msg=msg)
    if (istatus /= STATUS_INPUT_VALID) goto 100

    call curfit_get_wrk_size (m, k, nwrk, nest)

    if (size(knots) < nest .or. size(coefs) < nest) then
        if (present(msg)) msg = "knots and coefs arrays must be of length nest"
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
    else
        ! allocate real working array such that we can append
        ! weights at the end
        call ptr_work%assert_allocated (nrwrk=nwrk+m, niwrk=nest)
        ptr_work%rwrk(nwrk+1:nwrk+m) = 1.0_PREC
    end if

    ! call fitpack routine
    if (present(w)) then
        call fitpack_curfit_real64 (liopt, m, x, y, w, lxb, lxe, lk, &
            ls, nest, n, knots, coefs, lssr, &
            ptr_work%rwrk(1:nwrk), nwrk, ptr_work%iwrk, lstatus)
    else
        call fitpack_curfit_real64 (liopt, m, x, y, ptr_work%rwrk(nwrk+1:nwrk+m), &
            lxb, lxe, lk, ls, nest, n, knots, coefs, lssr, &
            ptr_work%rwrk(1:nwrk), nwrk, ptr_work%iwrk, lstatus)
    end if

    if (present(status)) status = lstatus
    if (present(ssr)) ssr = lssr

    return

100 if(present(status)) status = INTERP_STATUS_INVALID_INPUT

end subroutine

!-------------------------------------------------------------------------------
! CONCON wrapper routine

subroutine concon_wrapper (iopt, x, y, w, v, s, maxtr, maxbin, n, knots, coefs, &
        sx, bind, ssr, work, status)
    integer :: iopt, n, status, maxtr, maxbin
    real (PREC), dimension(:), contiguous :: x, y, w, v, knots, coefs, sx
    logical, dimension(:), contiguous :: bind
    real (PREC) :: s, ssr
    class (workspace), intent(in out), optional, target :: work

    intent(in) :: iopt, x, y, w, s, maxtr, maxbin
    intent(in out) :: knots, coefs, sx, bind, n, status, ssr
    optional :: iopt, w, s, maxtr, maxbin, bind, sx, status, ssr
    target :: w, sx, bind

    real (PREC), dimension(:), pointer, contiguous :: ptr_w => null(), ptr_sx => null()
    logical, dimension(:), pointer, contiguous :: ptr_bind => null()
    integer :: m, liopt, lstatus, nest, istatus, lmaxtr, lmaxbin
    integer :: nrwrk, rwrk_offset, niwrk, nlwrk, nrwrk_tot
    character (len=100) :: msg
    real (PREC) :: ls, lssr
    type (workspace), target :: lwork
    class (workspace), pointer :: ptr_work

    nullify (ptr_work)
    ! check for conformable arrays
    lstatus = INTERP_STATUS_INVALID_INPUT

    ! default values
    liopt = 0
    m = size(x)
    ls = 0.0_PREC
    ! lower bound of recommended interval for s
    if (present(w)) ls = m - sqrt(2*real(m, PREC))

    ! override default values if any of the optional parameters were provided
    if (present(iopt)) liopt = iopt
    if (present(s)) ls = s

    call check_input (x, y, w, knots=knots, coefs=coefs, v=v, sx=sx, &
        bind=bind, status=istatus, msg=msg)
    if (istatus /= STATUS_INPUT_VALID) goto 100

    nest = size(knots)
    lmaxbin = nest - 6
    lmaxtr = 200
    if (present(maxbin)) lmaxbin = maxbin
    if (present(maxtr)) lmaxtr = maxtr

    call concon_get_wrk_size (m, lmaxtr, lmaxbin, nest, nrwrk, niwrk)

    if (present(work)) then
        ptr_work => work
    else
        ptr_work => lwork
    end if

    ! work size for logical array
    nlwrk = 0
    ! offset for real work array; initially set to end of work array not
    ! including slices used for weights and sx
    rwrk_offset = nrwrk + 1
    nrwrk_tot = nrwrk
    if (.not. present(bind)) nlwrk = nlwrk + nest
    ! need to allocate workspace for uniform weights
    if (.not. present(w)) nrwrk_tot = nrwrk_tot + m
    ! need to allocate workspace for spline eval
    if (.not. present(sx)) nrwrk_tot = nrwrk_tot + m

    call ptr_work%assert_allocated (nrwrk=nrwrk_tot, niwrk=niwrk, nlwrk=nlwrk)

    ! assert that working arrays are sufficiently large
    if (present(w)) then
        ptr_w => w
    else
        ! allocate real working array such that we can append uniform
        ! weights at the end
        ptr_w => ptr_work%rwrk(rwrk_offset:rwrk_offset + m - 1)
        ptr_w = 1.0_PREC
        rwrk_offset = rwrk_offset + m
    end if

    if (present(sx)) then
        ptr_sx => sx
    else
        ! use sx allocated in working array that can be discarded after
        ! function call
        ptr_sx => ptr_work%rwrk(rwrk_offset:rwrk_offset + m - 1)
        rwrk_offset = rwrk_offset + m
    end if

    if (present(bind)) then
        ptr_bind => bind
    else
        ptr_bind => ptr_work%lwrk(1:nlwrk)
    end if

    ! initialize to zero as this is going to be used to store a tree where
    ! 0 denotes no node
    if (liopt == 0) ptr_work%iwrk = 0

    ! call fitpack routine wrapper
    call fitpack_concon_real64 (liopt, m, x, y, ptr_w, v, ls, nest, lmaxtr, lmaxbin, n, &
        knots, coefs, lssr, ptr_sx, ptr_bind, ptr_work%rwrk(1:nrwrk), nrwrk, &
        ptr_work%iwrk, niwrk, lstatus)

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

    call fitpack_splev_real64 (knots, n, coefs, lk, x, y, m, lext, lstatus)

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

    lk = DEFAULT_SPLINE_DEGREE
    lorder = 1
    lext = INTERP_EVAL_EXTRAPOLATE
    n = size(knots)
    m = size(x)

    if (present(k)) lk = k

    call check_input (x, y, k=lk, knots=knots, coefs=coefs, status=istatus, msg=msg)
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

    call fitpack_splder_real64 (knots, n, coefs, lk, lorder, x, y, m, lext, ptr_work%rwrk, lstatus)

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

pure function curfit_get_nest (m, k) result(res)
    integer, intent(in) :: k, m
    integer :: res

    if (k < MIN_SPLINE_DEGREE .or. k > MAX_SPLINE_DEGREE) then
        res = -1
    else
        res = max(m + k + 1, 2*k + 3)
    end if
end function

pure subroutine curfit_get_wrk_size (m, k, nwrk, nest)
    integer, intent(in) :: m, k
    integer, intent(out) :: nwrk, nest

    nest = curfit_get_nest (m, k)
    nwrk = m * (k + 1) + nest * (7 + 3 * k)
end subroutine

pure subroutine concon_get_wrk_size (m, maxtr, maxbin, nest, nrwrk, niwrk)
    integer :: m, maxtr, maxbin, nest, nrwrk, niwrk

    intent (in) :: m, maxtr, maxbin, nest
    intent (out) :: nrwrk, niwrk

    nrwrk = m * 4 + nest * 8 + maxbin * (maxbin + nest + 1)
    niwrk = maxtr * 4 + 2 * (maxbin + 1)
end subroutine

pure function concon_get_nest (m) result(nest)
    integer, intent(in) :: m
    integer :: nest

    ! return maximally required array size
    nest = m + 4
end function

end module
