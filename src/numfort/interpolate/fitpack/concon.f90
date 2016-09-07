module numfort_interpolate_concon

    use iso_fortran_env
    use numfort_core, only: comb
    use numfort_common, only: workspace, ENUM_KIND
    use numfort_interpolate_common
    use numfort_fitpack_interfaces
    use numfort_interpolate_result

    implicit none
    private

    integer, parameter :: PREC = real64

    ! used for internal input checking
    integer, parameter :: STATUS_INPUT_VALID = 0

    interface concon
        module procedure concon_wrapper
    end interface

    public :: concon_get_nest, concon

contains

! ******************************************************************************
! INPUT CHECKING used by various routines.

pure subroutine check_input (x, y, w, v, knots, coefs, sx, &
        bind, stat, msg)

    real (PREC), dimension(:) :: x, y, w, v, knots, coefs, sx
    logical, dimension(:) :: bind
    integer :: stat
    character (len=*) :: msg

    intent (in) :: x, y, w, v, knots, coefs, sx, bind
    intent (out) :: stat, msg
    optional :: w, sx, bind

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

    if (size(x) /= size(v)) then
        msg = "Non-conformable arrays x, v"
        goto 100
    end if

    if (size(knots) /= size(coefs)) then
        msg = "knots and coefs arrays must be of equal size"
        goto 100
    end if

    ! check for number of knots even though this is checked in concon again,
    ! as this defines nest, and nest is needed to compute maxtr
    if (size(knots) < 8) then
        msg = "knots/coefs arrays too small: size(knots) >= 8 required"
        goto 100
    end if

    if (present(sx)) then
        if (size(sx) /= size(x)) then
            msg = "x and sx must be of equal size"
            goto 100
        end if
    end if

    if (present(bind)) then
        if (size(bind) /= size(knots)) then
            msg = "knots and bind arrays must be of equal size"
            goto 100
        end if
    end if

    return

100 stat = INTERP_STATUS_INVALID_INPUT

end subroutine

! ******************************************************************************
! CONCON fitting procedure

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

    procedure (concon_if) :: concon

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

    call check_input (x, y, w, v, knots, coefs, sx, bind, istatus, msg)
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
    call concon (liopt, m, x, y, ptr_w, v, ls, nest, lmaxtr, lmaxbin, n, &
        knots, coefs, lssr, ptr_sx, ptr_bind, ptr_work%rwrk(1:nrwrk), nrwrk, &
        ptr_work%iwrk, niwrk, lstatus)

    if (present(status)) status = lstatus
    if (present(ssr)) ssr = lssr

    return

100 if(present(status)) status = INTERP_STATUS_INVALID_INPUT
    write (ERROR_UNIT, *) msg

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
