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
        module procedure splev_real64, splev_scalar_real64
    end interface

    interface curfit
        module procedure curfit_real64
    end interface

    interface concon
        module procedure concon_real64
    end interface

    interface splder
        module procedure splder_real64, splder_scalar_real64
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
    real (PREC), intent(in), dimension(:) :: x, y, knots, coefs
    integer, intent(in), optional :: k
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

pure subroutine check_eval_input (knots, coefs, n, k, x, y, order, ext, status, msg)
    real (PREC), intent(in), dimension(:) :: x, y, knots, coefs
    integer, intent(in), optional :: n
    integer, intent(in), optional :: k
    integer, intent(in), optional :: order
    integer, intent(in), optional :: ext

    integer (ENUM_KIND), intent(out) :: status
    character (len=*), intent(out), optional :: msg

    integer :: lk

    status = STATUS_INPUT_VALID

    lk = DEFAULT_SPLINE_DEGREE
    if (present(k)) lk = k

    if (size(knots) /= size(coefs)) then
        if (present(msg)) msg = 'knots and coefs array lengths differ'
        goto 100
    end if

    if (present(n)) then
        if (size(knots) < n) then
            if (present(msg)) &
                msg = "Length n too large for given knots/coefs arrays"
            goto 100
        end if
    end if

    ! check for supported spline degree
    if (lk < MIN_SPLINE_DEGREE .or. lk > MAX_SPLINE_DEGREE) then
        if (present(msg)) msg = "Unsupported spline degree k"
        goto 100
    end if

    if (size(x) /= size(y)) then
        if (present(msg)) msg = "x and y arrays must be of same length"
        goto 100
    end if

    if (present(order)) then
        if (order < 0 .or. order > lk) then
            if (present(msg)) msg = "Derivative order must satisfy 0<=order<=k"
        end if
    end if

    if (present(ext)) then
        call check_input_ext (ext, status, msg)
        if (status /= STATUS_INPUT_VALID) goto 100
    end if

    return

100 continue
    status = INTERP_STATUS_INVALID_INPUT

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

pure subroutine curfit_real64 (x, y, k, s, n, knots, coefs, &
        iopt, w, xb, xe, work, ssr, status, msg)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: y
    integer, intent(in), optional :: k
    real (PREC), intent(in), optional :: s
    integer, intent(out) :: n
    real (PREC), intent(in out), dimension(:), contiguous :: knots
    real (PREC), intent(out), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: iopt
    real (PREC), intent(in), dimension(:), contiguous, optional :: w
    real (PREC), intent(in), optional :: xe
    real (PREC), intent(in), optional :: xb
    class (workspace), intent(in out), optional, target :: work
    real (PREC), intent(out), optional :: ssr
    integer (ENUM_KIND), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: m, liopt, lk, nwrk, nest, lstatus
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

    call check_input (x, y, w, lk, knots, coefs, status=lstatus, msg=msg)
    if (lstatus /= STATUS_INPUT_VALID) goto 100

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

pure subroutine concon_real64 (x, y, v, s, n, &
        knots, coefs, iopt, w, maxtr, maxbin, work, ssr, status, msg)

    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: y
    real (PREC), intent(in), dimension(:), contiguous :: v
    real (PREC), intent(in), optional :: s
    integer, intent(out) :: n
    real (PREC), intent(in out), dimension(:), contiguous :: knots
    real (PREC), intent(out), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: iopt
    real (PREC), intent(in), dimension(:), contiguous, optional :: w
    integer, intent(in), optional :: maxtr
    integer, intent(in), optional :: maxbin
    class (workspace), intent(in out), optional, target :: work
    real (PREC), intent(out), optional :: ssr
    integer (ENUM_KIND), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: m, liopt, nest, lstatus, lmaxtr, lmaxbin
    integer :: nrwrk, jrwrk, niwrk, nlwrk, nrwrk_tot
    real (PREC) :: ls, lssr
    type (workspace), target :: lwork
    class (workspace), pointer :: ptr_work
    real (PREC), dimension(:), pointer, contiguous :: ptr_sx, ptr_w
    logical, dimension(:), pointer, contiguous :: ptr_bind

    nullify (ptr_work, ptr_sx, ptr_w, ptr_bind)

    call check_input (x, y, w, knots=knots, coefs=coefs, v=v, status=lstatus, msg=msg)
    if (lstatus /= STATUS_INPUT_VALID) goto 100

    ! check for conformable arrays
    lstatus = INTERP_STATUS_INVALID_INPUT

    ! default values
    liopt = 0
    m = size(x)
    ls = 0.0_PREC
    ! override default values if any of the optional parameters were provided
    ! lower bound of recommended interval for s
    if (present(w)) ls = m - sqrt(2*real(m, PREC))
    if (present(s)) ls = s
    if (present(iopt)) liopt = iopt

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
    nrwrk_tot = nrwrk
    ! elements containing bind array
    nlwrk = nlwrk + nest
    ! need to allocate workspace for uniform weights
    if (.not. present(w)) nrwrk_tot = nrwrk_tot + m
    ! need to allocate workspace for sx
    nrwrk_tot = nrwrk_tot + m

    call ptr_work%assert_allocated (nrwrk=nrwrk_tot, niwrk=niwrk, nlwrk=nlwrk)

    ! allocate default weights and sx after the actual real working array
    jrwrk = nrwrk + 1
    ! assert that working arrays are sufficiently large
    if (.not. present(w)) then
        ! allocate real working array such that we can append uniform
        ! weights at the end
        ptr_w => ptr_work%rwrk(jrwrk:jrwrk + m - 1)
        ptr_w = 1.0_PREC
        jrwrk = jrwrk + m
    end if

    ptr_sx => ptr_work%rwrk(jrwrk:jrwrk + m - 1)
    ptr_bind => ptr_work%lwrk(1:nest)

    ! initialize to zero as this is going to be used to store a tree where
    ! 0 denotes no node
    if (liopt == 0) ptr_work%iwrk = 0

    ! call fitpack routine wrapper
    if (present(w)) then
        call fitpack_concon_real64 (liopt, m, x, y, w, v, ls, nest, lmaxtr, &
            lmaxbin, n, knots, coefs, lssr, ptr_sx, ptr_bind, &
            ptr_work%rwrk(1:nrwrk), nrwrk, ptr_work%iwrk, niwrk, lstatus)
    else
        call fitpack_concon_real64 (liopt, m, x, y, ptr_w, v, ls, nest, lmaxtr, &
            lmaxbin, n, knots, coefs, lssr, ptr_sx, ptr_bind, &
            ptr_work%rwrk(1:nrwrk), nrwrk, ptr_work%iwrk, niwrk, lstatus)
    end if

    if (present(status)) status = lstatus
    if (present(ssr)) ssr = lssr

    return

100 continue
    if(present(status)) status = INTERP_STATUS_INVALID_INPUT

end subroutine

! ******************************************************************************
! CURFIT evaluation

pure subroutine splev_real64 (knots, coefs, n, k, x, y, ext, status)
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: n
    integer, intent(in), optional :: k
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer, intent(in), optional :: ext
    integer, intent(out), optional :: status

    integer :: m, ln, lstatus, lext, lk

    lstatus = INTERP_STATUS_INVALID_INPUT

    call check_eval_input (knots, coefs, n, k, x, y, ext=ext, status=lstatus)
    if (lstatus /= STATUS_INPUT_VALID) goto 100

    m = size(x)
    ! by default assume that all elements in knots/coefs array are valid
    ln = size(knots)
    ! assume cubic splines by default
    lk = DEFAULT_SPLINE_DEGREE
    ! by default set function values outside of domain to boundary values
    lext = INTERP_EVAL_EXTRAPOLATE

    ! If n is present, cap number of knots/coefs
    if (present(n)) ln = n
    if (present(k)) lk = k

    call fitpack_splev_real64 (knots, ln, coefs, lk, x, y, m, lext, lstatus)

    if (present(status)) status = lstatus

    return

100 continue
    if (present(status)) status = INTERP_STATUS_INVALID_INPUT

end subroutine

pure subroutine splev_scalar_real64 (knots, coefs, n, k, x, y, ext, status)
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: n
    integer, intent(in), optional :: k
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in), optional :: ext
    integer, intent(out), optional :: status

    real (PREC) :: lx(1), ly(1)

    lx(1) = x

    call splev_real64 (knots, coefs, n, k, lx, ly, ext, status)

    ! write back result
    y = ly(1)
end subroutine

! ******************************************************************************
! Evaluation of derivatives

pure subroutine splder_real64 (knots, coefs, n, k, order, x, y, ext, work, status)
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: n
    integer, intent(in), optional :: k
    integer, intent(in), optional :: order
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer, intent(in), optional :: ext
    class (workspace), intent(in out), optional, target :: work
    integer, intent(out), optional :: status

    type (workspace), target :: lwork
    type (workspace), pointer :: ptr_work

    integer :: lext, lorder, lstatus, lk, ln, m

    call check_eval_input (knots, coefs, n, k, x, y, order, ext, status=lstatus)
    if (lstatus /= STATUS_INPUT_VALID) goto 100

    lk = DEFAULT_SPLINE_DEGREE
    lorder = 1
    lext = INTERP_EVAL_EXTRAPOLATE
    ln = size(knots)
    m = size(x)

    if (present(n)) ln = n
    if (present(k)) lk = k

    if (present(work)) then
        ptr_work => work
    else
        ptr_work => lwork
    end if

    call ptr_work%assert_allocated (nrwrk=ln)

    call fitpack_splder_real64 (knots, ln, coefs, lk, lorder, x, y, m, &
        lext, ptr_work%rwrk, lstatus)

    if (present(status)) status = lstatus

    return

100 continue
    if (present(status)) status = INTERP_STATUS_INVALID_INPUT
end subroutine

! splder_scalar_real64 evaluates the derivate of a spline at a single scalar point
! x and stores the result in scalar y.
! Implemented as a wrapper for SPLDER_WRAPPER that creates local size-1 arrays.
pure subroutine splder_scalar_real64 (knots, coefs, n, k, order, x, y, ext, work, status)
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: n
    integer, intent(in), optional :: k
    integer, intent(in), optional :: order
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in), optional :: ext
    class (workspace), intent(in out), optional, target :: work
    integer, intent(out), optional :: status

    real (PREC), dimension(1) :: lx, ly

    lx(1) = x
    ly(1) = y

    call splder_real64 (knots, coefs, n, k, order, lx, ly, ext, work, status)

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
