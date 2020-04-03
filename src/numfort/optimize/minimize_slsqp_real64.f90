


module numfort_optimize_slsqp_real64

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_optimize_result_real64
    use numfort_optimize_interfaces_real64
    use numfort_optimize_fwrapper_real64
    use slsqp_mod, only: slsqp_orig => slsqp, slsqp_data => slsqp_data_real64, &
        linmin_data => linmin_data_real64

    implicit none

    integer, parameter :: PREC = real64

    private

    public :: minimize_slsqp

    interface minimize_slsqp
        module procedure minimize_slsqp, minimize_slsqp_args
    end interface

    integer, parameter :: MODE_INIT = 0
    integer, parameter :: MODE_EVAL_FUNCS = 1
    integer, parameter :: MODE_EVAL_JAC = -1

    character (*), parameter :: SLSQP_MSG_MAP(0:9) = [ &
            "Optimization terminated successfully               ", &
            "                                                   ", &
            "Move equality constraints than indepedent variables", &
            "More than 3*n iterations in LSQ subproblem         ", &
            "Incompatible inequality constraints                ", &
            "Singular matrix E in LSQ subproblem                ", &
            "Singular matrix C in LSQ subproblem                ", &
            "Rank-deficient equality constraint in HFTI         ", &
            "Positive directional derivative in line search     ", &
            "Iteration limit exceeded                           "  &
        ]
        !*  Mapping of SLSQP status codes to status messages. Note: no msg
        !   for value 1 as this status code is not used on termination.

    integer (NF_ENUM_KIND), parameter :: SLSQP_STATUS_MAP(0:9) = [ &
            NF_STATUS_OK, &
            NF_STATUS_UNDEFINED, &
            NF_STATUS_INVALID_ARG, &
            NF_STATUS_NOT_CONVERGED, &
            NF_STATUS_INVALID_STATE, &
            NF_STATUS_INVALID_STATE, &
            NF_STATUS_INVALID_STATE, &
            NF_STATUS_INVALID_STATE, &
            NF_STATUS_INVALID_STATE, &
            NF_STATUS_MAX_ITER &
        ]
        !*  Mapping of SLSQP to numfort status codes. Note: value 1
        !   is never returned on termination.

    contains



subroutine minimize_slsqp (fobj, x, lbounds, ubounds, m, meq, f_eqcons, &
        f_ieqcons, maxiter, linesearch, tol, work, res)

    procedure (fvs_fcn_jac_opt) :: fobj
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: m, meq
    procedure (fvv_fcn_jac_opt), optional :: f_eqcons
    procedure (fvv_fcn_jac_opt), optional :: f_ieqcons
    integer, intent(in), optional :: maxiter
    integer (NF_ENUM_KIND), intent(in), optional :: linesearch
    real (PREC), intent(in), optional :: tol
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vs) :: fobj_wrapper
    type (fwrapper_vv) :: f_eqcons_wrapper, f_ieqcons_wrapper

    call wrap_procedure (fobj_wrapper, fcn_jac_opt=fobj)
    call wrap_procedure (f_eqcons_wrapper, fcn_jac_opt=f_eqcons)
    call wrap_procedure (f_ieqcons_wrapper, fcn_jac_opt=f_ieqcons)

    call minimize_slsqp_impl (fobj_wrapper, x, lbounds, ubounds, m, meq, &
        f_eqcons_wrapper, f_ieqcons_wrapper, maxiter, linesearch, tol, work, res)

end subroutine


subroutine minimize_slsqp_args (fobj, x, args, lbounds, ubounds, m, meq, f_eqcons, &
        f_ieqcons, maxiter, linesearch, tol, work, res)

    procedure (fvs_fcn_jac_opt_args) :: fobj
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: m, meq
    procedure (fvv_fcn_jac_opt_args), optional :: f_eqcons
    procedure (fvv_fcn_jac_opt_args), optional :: f_ieqcons
    integer, intent(in), optional :: maxiter
    integer (NF_ENUM_KIND), intent(in), optional :: linesearch
    real (PREC), intent(in), optional :: tol
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vs) :: fobj_wrapper
    type (fwrapper_vv) :: f_eqcons_wrapper, f_ieqcons_wrapper

    call wrap_procedure (fobj_wrapper, fcn_jac_opt_args=fobj, args=args)
    call wrap_procedure (f_eqcons_wrapper, fcn_jac_opt_args=f_eqcons, args=args)
    call wrap_procedure (f_ieqcons_wrapper, fcn_jac_opt_args=f_ieqcons, args=args)

    call minimize_slsqp_impl (fobj_wrapper, x, lbounds, ubounds, m, meq, &
        f_eqcons_wrapper, f_ieqcons_wrapper, maxiter, linesearch, tol, work, res)

end subroutine



subroutine minimize_slsqp_impl (fobj, x, lbounds, ubounds, m, meq, f_eqcons, &
        f_ieqcons, maxiter, linesearch, tol, work, res)

    type (fwrapper_vs) :: fobj
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: m, meq
    type (fwrapper_vv) :: f_eqcons
    type (fwrapper_vv) :: f_ieqcons
    integer, intent(in), optional :: maxiter
    integer (NF_ENUM_KIND), intent(in), optional :: linesearch
    real (PREC), intent(in), optional :: tol
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (optim_result), pointer :: ptr_res
        !   Pointer to local or user-provided result object

    integer :: iter, lmaxiter
    integer (NF_ENUM_KIND) :: llinesearch
    real (PREC) :: acc

    type (slsqp_data) :: dat
        !   Stores persistent internal data used in SLSQPB routine
    type (linmin_data) :: dat_lm
        !   Stores persistent internal data used in LINMIN routine
    type (workspace), pointer :: ptr_work
    real (PREC), dimension(:), pointer, contiguous :: ptr_xlb, ptr_xub
    real (PREC), dimension(:), pointer, contiguous :: ptr_g, ptr_w, ptr_c
    real (PREC), dimension(:), pointer, contiguous :: ptr_cx_eq, ptr_cx_ieq
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_cpx_eq, ptr_cpx_ieq
        !   Pointers to arrays that store Jacobians of eq. and ineq. constraints
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_a
        !   Pointer to array that stores stacked Jacobians of eq. and ineq.
        !   constraints.
    integer :: nrwrk, niwrk, ioffset
    integer, dimension(2) :: shp

    integer :: n, n1, mineq, mode, lmeq, lm, mieq, k, lda, lw
    logical :: has_eq, has_ieq
    real (PREC) :: fx, ltol
    type (status_t) :: status

    status = NF_STATUS_OK

    nullify (ptr_work)
    nullify (ptr_xlb, ptr_xub)
    nullify (ptr_g, ptr_c, ptr_a, ptr_w)
    nullify (ptr_cx_eq, ptr_cx_ieq)
    nullify (ptr_cpx_eq, ptr_cpx_ieq)

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call slsqp_check_input (x, lbounds, ubounds, m ,meq, f_eqcons, f_ieqcons, &
        maxiter, linesearch, tol, status, ptr_res%msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    ! Determine if enough information is present for equality or inequality
    ! constraints.
    has_eq = is_associated(f_eqcons) .and. (present(meq) .or. &
        (.not. is_associated(f_ieqcons) .and. present(m)))

    has_ieq = (is_associated(f_ieqcons) .and. present(m)) .and. &
        (.not. is_associated(f_eqcons) .or. present(meq))

    lmaxiter = 100
    ltol = 1.0d-6
    lmeq = 0
    lm = 0
    llinesearch = NF_LINESEARCH_BACKTRACK
    if (present(maxiter)) lmaxiter = maxiter
    if (present(tol)) ltol = tol
    if (present(m)) lm = m
    if (present(meq)) lmeq = meq
    if (present(linesearch)) llinesearch = linesearch

    ! Default initialize some values that might be retured as attributes
    ! of OPTIM_RESULT before minimization starts if an error is encountered
    fx = 0.0
    iter = 0

    ! if equality constraint specified, but no inequality constraint,
    ! m = meq and allow for missing meq if m is present.
    if (has_eq .and. .not. is_associated(f_ieqcons)) then
        if (present(meq) .and. .not. present(m)) then
            lm = lmeq
        else if (present(m) .and. .not. present(meq)) then
            lmeq = lm
        end if
    end if

    ! if inequality constraint specified, but no equality constraint,
    ! meq = 0 and mieq = m, thus allow for missing meq
    if (has_ieq .and. .not. is_associated(f_eqcons)) lmeq = 0

    ! implied number of inequality constraints
    mieq = lm - lmeq

    n = size(x)
    mineq = mieq + 2 * (n + 1)

    niwrk = mineq
    n1 = n + 1
    ! Workspace used directly by SLSQP: Take calculations from implementation.
    lw = (3*n1+lm)*(n1+1)+(n1-lmeq+1)*(mineq+2) + 2*mineq+(n1+mineq)*(n1-lmeq) &
        + 2*lmeq + n1 + ((n+1)*n)/2 + 2*lm + 3*n + 3*n1 + 1
    ! potential bug in underlying SLSQP routine that leads to random memory
    ! errors. Maybe WA is accessed past its boundary, so increase by 1 KiB for
    ! now and hope for the best.
    lw = lw + 128
    nrwrk = lw

    ! Add space for lower / upper bounds
    nrwrk = nrwrk + 2 * n
    ! Add space for gradient (SLSQP expects this to be size n+1 for whatever reason)
    nrwrk = nrwrk + n + 1
    ! Add space for eq. and ineq. constraint values; SLSQP expects constraint
    ! array to be at least size 1 regardless of how many constraints are present.
    nrwrk = nrwrk + max(1, lm)
    ! Space for stacked constraint Jacobian matrix A. These are just the
    ! "vertically" stacked Jacobians of eq. and ineq. constraints, but
    ! due to vertical stacking we need separate arrays to pass to
    ! f_eqcons and f_ieqcons if we want arguments to be contiguous. Hence
    ! allocate additional set of arrays.
    nrwrk = nrwrk + max(1, lm) * (n+1)
    ! Space for Jacobians of eq. and ineq. constraints
    nrwrk = nrwrk + lmeq * n + mieq * n

    ! Allocate workspace arrays
    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal state in workspace object, in particular index offsets
    ! (this does not deallocate working arrays)
    call workspace_reset (ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk, niwrk=niwrk)

    call workspace_get_ptr (ptr_work, n, ptr_xlb, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, ptr_xub, status)
    if (status /= NF_STATUS_OK) goto 100

    ! Assert that bounds values are as expected by underlying routine
    call slsqp_sanitize_bounds (x, lbounds, ubounds, ptr_xlb, ptr_xub, ptr_res)
    if (NF_STATUS_INVALID_ARG .in. ptr_res%status) goto 100
    ! Clip initial guess such that it satisfied any lower/upper bounds
    call slsqp_clip_init (x, ptr_xlb, ptr_xub)

    ! initial value for mode parameter
    mode = 0

    ! SLSQP implementation expects gradient array to be of size (n+1)
    call workspace_get_ptr (ptr_work, n+1, ptr_g, status)
    if (status /= NF_STATUS_OK) goto 100

    ! Pointers to data containing evaluated constraints
    if (.not. has_eq .and. .not. has_ieq) then
        ! Array C is expected to be at least size 1
        call workspace_get_ptr (ptr_work, 1, ptr_c, status)
        if (status /= NF_STATUS_OK) goto 100

        ! Add "empty" row of stacked Jacobian matrices as SLSQP expects
        ! argument A to be at least (1, n+1).
        ! Assign shape such that gfortran does not create warnings about
        ! temporaray arrays.
        shp(1) = 1
        shp(2) = n+1
        call workspace_get_ptr (ptr_work, shp, ptr_a, status)
        if (status /= NF_STATUS_OK) goto 100
    else
        ! Array C is expected to contain eq. and ineq. constraint values
        ! concatenated together.
        call workspace_get_ptr (ptr_work, lm, ptr_c, status)
        if (status /= NF_STATUS_OK) goto 100

        shp = lm
        shp(2) = n+1
        call workspace_get_ptr (ptr_work, shp, ptr_a, status)
        if (status /= NF_STATUS_OK) goto 100

        ioffset = 0
        if (has_eq) then
            ptr_cx_eq => ptr_c(1:lmeq)
            ioffset = lmeq

            shp(1) = lmeq
            shp(2) = n
            call workspace_get_ptr (ptr_work, shp, ptr_cpx_eq, status)
            if (status /= NF_STATUS_OK) goto 100
        end if
        if (has_ieq) then
            ptr_cx_ieq(1:mieq) => ptr_c(ioffset+1:lm)

            shp(1) = mieq
            shp(2) = n
            call workspace_get_ptr (ptr_work, shp, ptr_cpx_ieq, status)
            if (status /= NF_STATUS_OK) goto 100
        end if
    end if

    ! Assign remaining arguments
    lda = max(1, lm)
    ! Pointer to workspace used directly by SLSQP
    call workspace_get_ptr (ptr_work, lw, ptr_w, status)
    if (status /= NF_STATUS_OK) goto 100

    ! ACC argument: for exact linesearch this has to be the negative tolerance
    acc = ltol
    if (llinesearch == NF_LINESEARCH_EXACT) acc = - abs(ltol)

    ! Overwrite any garbage. In particular, column n+1 of these arrays
    ! will never be used to store data, so no idea what SLSQP is doing with in.
    ptr_g(:) = 0.0_PREC
    ptr_a(:,:) = 0.0_PREC

    ! Initial value passed into SLSQP determines max. iterations, but
    ! value is subsequently reset and incremented as new iteration starts.
    ! SLSQP aborts by itself with the appropriate exit code when max. iter.
    ! count is exceeded.
    iter = lmaxiter

    do while (.true.)
        if (mode == MODE_INIT .or. mode == MODE_EVAL_FUNCS) then
            ! Compute objective function fx
            call dispatch (fobj, x, fx)

            ! Compute equality constraints
            if (has_eq) then
                call dispatch (f_eqcons, x, ptr_cx_eq)
            end if

            ! Compute inequality constraints
            if (has_ieq) then
                call dispatch (f_ieqcons, x, ptr_cx_ieq)
            end if
        end if

        if (mode == MODE_INIT .or. mode == MODE_EVAL_JAC) then
            ! Compute gradient of objective function
            call dispatch_jac (fobj, x, fpx=ptr_g(1:n))

            ioffset = 0
            if (has_eq) then
                ! Compute n-by-meq Jacobian of equality constraints
                call dispatch_jac (f_eqcons, x, fpx=ptr_cpx_eq)
                ! Copy Jacobian into array A using a loop, otherwise
                ! (at least) gfortran allocates a temporary array.
                forall (k=1:lmeq) ptr_a(k,1:n) = ptr_cpx_eq(k,:)
                ioffset = lmeq
            end if

            if (has_ieq) then
                call dispatch_jac (f_ieqcons, x, fpx=ptr_cpx_ieq)
                ! Concatenate Jacobian into array A, taking into account
                ! whether Jacobian of eq. constr. is already present
                forall (k=1:mieq) ptr_a(ioffset+k,1:n) = ptr_cpx_ieq(k,:)
            end if
        end if

        ! call original SLSQP routine
        call slsqp_orig (dat, dat_lm, lm, lmeq, lda, n, x, ptr_xlb, ptr_xub, fx, &
            ptr_c, ptr_g, ptr_a, acc, iter, mode, ptr_w, lw, &
            ptr_work%iwrk, mineq)

        ! Exit in all cases when there is an error or convergence was achieved.
        if (mode /= MODE_EVAL_FUNCS .and. mode /= MODE_EVAL_JAC) then
            ! Update status code and status message on termination
            status = SLSQP_STATUS_MAP(mode)
            status%code_orig = mode
            call result_update (ptr_res, msg=SLSQP_MSG_MAP(mode))
            goto 100
        end if
    end do


100 continue

    ! Update result object only if we are returning to client code
    if (present(res)) then
        call result_update (ptr_res, x, fx, nit=iter, nfev=fobj%nfev, &
            status=status)
    end if

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

end subroutine



subroutine slsqp_check_input (x, lbounds, ubounds, m, meq, &
        f_eqcons, f_ieqcons, maxiter, linesearch, tol, status, msg)

    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds, ubounds
    integer, intent(in), optional :: m, meq
    type (fwrapper_vv), intent(in) :: f_eqcons
    type (fwrapper_vv), optional :: f_ieqcons
    integer, intent(in), optional :: maxiter
    integer (NF_ENUM_KIND), intent(in), optional :: linesearch
    real (PREC), intent(in), optional :: tol
    type (status_t), intent(inout) :: status
    character (*), intent(out) :: msg

    integer, parameter :: LS_VALID(2) = [NF_LINESEARCH_BACKTRACK, NF_LINESEARCH_EXACT]
        !   Valid values for linesearch argument
    logical :: has_eq, has_ieq
    integer :: n

    status = NF_STATUS_OK

    n = size(x)

    if (present(lbounds)) then
        if (size(lbounds) < n) then
            msg = "Argument 'lbounds': non-conformable array size"
            goto 100
        end if
    end if

    if (present(ubounds)) then
        if (size(ubounds) < n) then
            msg = "Argument 'ubounds': non-conformable array size"
            goto 100
        end if
    end if

    has_eq = is_associated(f_eqcons) .and. (present(meq) .or. &
        (.not. is_associated(f_ieqcons) .and. present(m)))

    has_ieq = (is_associated(f_ieqcons) .and. present(m)) .and. &
        (.not. is_associated(f_eqcons) .or. present(meq))

    if (is_associated(f_eqcons) .and. .not. has_eq) then
        msg = "Number of equality constraints not specified"
        goto 100
    end if

    if (is_associated(f_ieqcons) .and. .not. has_ieq) then
        msg = "Number of inequality constraints not specified"
        goto 100
    end if

    call check_positive (1, meq, 'meq', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    call check_positive (1, m, 'm', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    call check_positive (1, maxiter, 'maxiter', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    call check_enum (linesearch, LS_VALID, 'linesearch', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    call check_positive (1.0_PREC, tol, 'tol', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



subroutine slsqp_sanitize_bounds (x, xlb_in, xub_in, xlb, xub, res)

    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:), optional :: xlb_in, xub_in
    real (PREC), intent(out), dimension(:) :: xlb, xub
    type (optim_result), intent(inout) :: res

    integer :: n, i
    real (PREC) :: POS_INF, NEG_INF
    type (status_t) :: status

    n = size(x)

    POS_INF = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    NEG_INF = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)

    xlb = NEG_INF
    xub = POS_INF

    if (present(xlb_in)) then
        do i = 1, n
            if (ieee_is_finite (xlb_in(i))) then
                xlb(i) = xlb_in(i)
            else
                xlb(i) = NEG_INF
            end if
        end do
    end if

    if (present(xub_in)) then
        do i = 1, n
            if (ieee_is_finite (xub_in(i))) then
                xub(i) = xub_in(i)
            else
                xub(i) = POS_INF
            end if
        end do
    end if

    ! Check that xlb <= xub holds
    do i = 1, n
        if (xub(i) < xlb(i)) then
            status = NF_STATUS_INVALID_ARG
            call result_update (res, status=status, &
                msg="Invalid lower and/or upper bound specified")
            return
        end if
    end do

end subroutine



subroutine slsqp_clip_init (x, xlb, xub)

    real (PREC), intent(inout), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: xlb, xub

    integer :: i

    do i = 1, size(x)
        if (ieee_is_finite (xlb(i))) then
            x(i) = max(xlb(i), x(i))
        end if
        if (ieee_is_finite (xub(i))) then
            x(i) = min(xub(i), x(i))
        end if
    end do

end subroutine

end module
