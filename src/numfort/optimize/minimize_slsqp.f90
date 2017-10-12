module numfort_optimize_slsqp

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_common_workspace
    use numfort_optimize_result
    use slsqp_mod, slsqp_impl => slsqp

    implicit none

    private

    public :: slsqp

    interface slsqp
        module procedure slsqp_wrapper_real64
    end interface

    interface slsqp_check_input
        module procedure slsqp_check_input_real64
    end interface

    interface slsqp_clip_init
        module procedure slsqp_clip_init_real64
    end interface

    interface slsqp_sanitize_bounds
        module procedure slsqp_sanitize_bounds_real64
    end interface

    integer, parameter :: MODE_INIT = 0
    integer, parameter :: MODE_EVAL_FUNCS = 1
    integer, parameter :: MODE_EVAL_JAC = -1

    abstract interface
        subroutine fobj_real64 (x, fx, fpx, args)
            import real64
            real (real64), intent(in), dimension(:) :: x
            real (real64), intent(out), optional :: fx
            real (real64), intent(out), dimension(:), optional :: fpx
            real (real64), intent(in), dimension(:), optional :: args
        end subroutine

        subroutine f_eqcons_real64 (x, fx, fpx, args)
            import real64
            real (real64), intent(in), dimension(:) :: x
            real (real64), intent(out), dimension(:), optional :: fx
            real (real64), intent(out), dimension(:,:), optional :: fpx
            real (real64), intent(in), dimension(:), optional :: args
        end subroutine

        subroutine f_ieqcons_real64 (x, fx, fpx, args)
            import real64
            real (real64), intent(in), dimension(:) :: x
            real (real64), intent(out), dimension(:), optional :: fx
            real (real64), intent(out), dimension(:,:), optional :: fpx
            real (real64), intent(in), dimension(:), optional :: args
        end subroutine
    end interface

contains


subroutine slsqp_wrapper_real64 (func, x, lbounds, ubounds, f_eqcons, meq, &
        f_ieqcons, mieq, args, maxiter, tol, work, res)

    integer, parameter :: PREC = real64

    procedure (fobj_real64) :: func
    real (PREC), intent(in out), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    procedure (f_eqcons_real64), optional :: f_eqcons
    integer, intent(in), optional :: meq
    procedure (f_ieqcons_real64), optional :: f_ieqcons
    integer, intent(in), optional :: mieq
    real (PREC), intent(in), dimension(:), optional :: args
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: tol
    type (workspace_real64), intent(in out), optional :: work
    type (optim_result_real64), intent(in out), optional :: res

    type (status_t) :: status
    character (100) :: msg
    integer :: iter, nfeval, lmaxiter

    type (slsqp_data_real64) :: dat
        !   Stores internal data used in SLSQPB that was originally declared
        !   with a SAVE attribute.
    type (workspace_real64), pointer :: ptr_work
    real (PREC), dimension(:), pointer, contiguous :: ptr_xlb, ptr_xub
    real (PREC), dimension(:), pointer, contiguous :: ptr_g, ptr_w, ptr_c, ptr_tmp
    real (PREC), dimension(:), pointer, contiguous :: ptr_cx_eq, ptr_cx_ieq
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_cpx_eq, ptr_cpx_ieq
        !   Pointers to arrays that store Jacobians of eq. and ineq. constraints
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_a
        !   Pointer to array that stores stacked Jacobians of eq. and ineq.
        !   constraints.
    integer :: nrwrk, niwrk, ioffset

    integer :: n, n1, m, mineq, mode, lmeq, lmieq, k, lda, lw
    logical :: has_eq, has_ieq
    real (PREC) :: fx, ltol

    status = NF_STATUS_INVALID_ARG
    nullify (ptr_work)
    nullify (ptr_xlb, ptr_xub)
    nullify (ptr_g, ptr_c, ptr_a, ptr_w, ptr_tmp)
    nullify (ptr_cx_eq, ptr_cx_ieq)
    nullify (ptr_cpx_eq, ptr_cpx_ieq)

    call slsqp_check_input (f_eqcons, meq, f_ieqcons, mieq, maxiter, tol, status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    has_eq = present(f_eqcons) .and. present(meq)
    has_ieq = present(f_ieqcons) .and. present(mieq)

    lmaxiter = 100
    ltol = 1.0d-6
    lmeq = 0
    lmieq = 0
    if (present(maxiter)) lmaxiter = maxiter
    if (present(tol)) ltol = tol
    ! If any equality or inequality constraints were specified, store dimension
    ! of constraint function range.
    if (has_eq) lmeq = meq
    if (has_ieq) lmieq = mieq

    n = size(x)
    m = lmeq + lmieq
    mineq = lmieq + 2 * (n + 1)

    niwrk = mineq
    n1 = n + 1
    ! Workspace used directly by SLSQP
    lw = (3*n1+m)*(n1+1)+(n1-lmeq+1)*(mineq+2) + 2*mineq+(n1+mineq)*(n1-lmeq) &
        + 2*lmeq + n1 + ((n+1)*n)/2 + 2*m + 3*n + 3*n1 + 1
    nrwrk = lw

    ! Add space for lower / upper bounds
    nrwrk = nrwrk + 2 * n
    ! Add space for gradient (SLSQP expects this to be size n+1 for whatever reason)
    nrwrk = nrwrk + n + 1
    ! Add space for eq. and ineq. constraint values; SLSQP expects constraint
    ! array to be at least size 1 regardless of how many constraints are present.
    nrwrk = nrwrk + max(1, m)
    ! Space for stacked constraint Jacobian matrix A. These are just the
    ! "vertically" stacked Jacobians of eq. and ineq. constraints, but
    ! due to vertical stacking we need separate arrays to pass to
    ! f_eqcons and f_ieqcons if we want arguments to be contiguous. Hence
    ! allocate additional set of arrays.
    nrwrk = nrwrk + max(1, m) * (n+1)
    ! Space for Jacobians of eq. and ineq. constraints
    nrwrk = nrwrk + lmeq * n + lmieq * n

    ! Allocate workspace arrays
    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal state in workspace object, in particular index offsets
    ! (this does not deallocate working arrays)
    call workspace_reset (ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk, niwrk=niwrk)

    ptr_xlb => workspace_get_ptr (ptr_work, n)
    ptr_xub => workspace_get_ptr (ptr_work, n)

    ! Assert that bounds values are as expected by underlying routine
    call slsqp_sanitize_bounds (lbounds, ubounds, ptr_xlb, ptr_xub)
    ! Clip initial guess such that it satisfied any lower/upper bounds
    call slsqp_clip_init (x, ptr_xlb, ptr_xub)

    ! initial value for mode parameter
    mode = 0

    ! SLSQP implementation expects gradient array to be of size (n+1)
    ptr_g => workspace_get_ptr (ptr_work, n+1)

    ! Pointers to data containing evaluated constraints
    if (.not. has_eq .and. .not. has_ieq) then
        ! Array C is expected to be at least size 1
        ptr_c => workspace_get_ptr (ptr_work, 1)

        ! Add "empty" row of stacked Jacobian matrices as SLSQP expects
        ! argument A to be at least (1, n+1).
        ! BUG: gfortran-5.x segfaults when assigning a pointer to a 1d array
        ! return from a function to a 2d array. Use temporary 1d-array pointer
        ! as workaround.
        ptr_tmp => workspace_get_ptr (ptr_work, n+1)
        ptr_a(1:1,1:n+1) => ptr_tmp
        nullify (ptr_tmp)
    else
        ! Array C is expected to contain eq. and ineq. constraint values
        ! concatenated together.
        ptr_c => workspace_get_ptr (ptr_work, m)
        ptr_tmp => workspace_get_ptr (ptr_work, m * (n+1))
        ptr_a(1:m,1:n+1) => ptr_tmp
        nullify (ptr_tmp)

        ioffset = 0
        if (has_eq) then
            ptr_cx_eq => ptr_c(1:lmeq)
            ioffset = lmeq

            k = lmeq * (n+1)
            ptr_tmp => workspace_get_ptr (ptr_work, k)
            ptr_cpx_eq(1:lmeq,1:n+1) => ptr_tmp
            nullify (ptr_tmp)
        end if
        if (has_ieq) then
            ptr_cx_ieq(1:n) => ptr_c(ioffset+1:m)

            k = lmieq * (n+1)
            ptr_tmp => workspace_get_ptr (ptr_work, k)
            ptr_cpx_ieq(1:lmieq,1:n+1) => ptr_tmp
            nullify (ptr_tmp)
        end if
    end if


    do iter = 1, lmaxiter
        if (mode == MODE_INIT .or. mode == MODE_EVAL_FUNCS) then
            ! Compute objective function fx
            call func (x, fx, args=args)
            ! Increment function evaluation counter
            nfeval = nfeval + 1

            ! Compute equality constraints
            if (has_eq) then
                call f_eqcons (x, ptr_cx_eq, args=args)
            end if

            ! Compute inequality constraints
            if (has_ieq) then
                call f_ieqcons (x, ptr_cx_ieq, args=args)
            end if
        end if

        if (mode == MODE_INIT .or. mode == MODE_EVAL_JAC) then
            ! Compute gradient of objective function
            call func (x, fpx=ptr_g(1:n), args=args)

            ioffset = 0
            if (has_eq) then
                ! Compute n-by-meq Jacobian of equality constraints
                call f_eqcons (x, fpx=ptr_cpx_eq, args=args)
                ! Copy Jacobian into array A using a loop, otherwise
                ! (at least) gfortran allocates a temporary array.
                forall (k=1:lmeq) ptr_a(k,:) = ptr_cpx_eq(k,:)
            end if

            if (has_ieq) then
                call f_ieqcons (x, fpx=ptr_cpx_ieq, args=args)
                ! Concatenate Jacobian into array A, taking into account
                ! whether Jacobian of eq. constr. is already present
                forall (k=ioffset+1:ioffset+lmieq) ptr_a(k,:) = ptr_cpx_ieq(k,:)
            end if
        end if

        ! Assign remaining arguments
        lda = max(1, m)
        ! Pointer to workspace used directly by SLSQP
        ptr_w => workspace_get_ptr (ptr_work, lw)
        ! call original SLSQP routine
        call slsqp_impl (dat, m, lmeq, lda, n, x, ptr_xlb, ptr_xub, fx, &
            ptr_c, ptr_g, ptr_a, ltol, lmaxiter, mode, ptr_w, lw, &
            ptr_work%iwrk, mineq)

        select case (mode)
        case (0)
            status = NF_STATUS_OK
            msg = "Optimization terminated successfully"
        case (2)
            status = NF_STATUS_INVALID_ARG
            msg = "Move equality constraints than indepedent variables"
        case (3)
            status = NF_STATUS_NOT_CONVERGED
            msg = "More than 3*n iterations in LSQ subproblem"
        case (4)
            status = NF_STATUS_INVALID_STATE
            msg = "Incompatible inequality constraints"
        case (5)
            status = NF_STATUS_INVALID_STATE
            msg = "Singular matrix E in LSQ subproblem"
        case (6)
            status = NF_STATUS_INVALID_STATE
            msg = "Singular matrix C in LSQ subproblem"
        case (7)
            status = NF_STATUS_INVALID_STATE
            msg = "Rank-deficient equality constraint in HFTI"
        case (8)
            status = NF_STATUS_INVALID_STATE
            msg = "Positive directional derivative in line searc"
        case (9)
            status = NF_STATUS_MAX_ITER
            msg = "Iteration limit exceeded"
        end select

        ! Exit in all cases when there is an error or convergence was achieved.
        if (mode /= MODE_EVAL_FUNCS .and. mode /= MODE_EVAL_JAC) goto 100
    end do


100 continue

    if (present(res)) then
        call result_update (res, x, fx, status, iter, nfeval, msg)
    end if

    call assert_dealloc_ptr (work, ptr_work)

end subroutine



subroutine slsqp_check_input_real64 (f_eqcons, meq, f_ieqcons, mieq, maxiter, &
        tol, status, msg)

    integer, parameter :: PREC = real64

    procedure (f_eqcons_real64), optional :: f_eqcons
    integer, intent(in), optional :: meq
    procedure (f_ieqcons_real64), optional :: f_ieqcons
    integer, intent(in), optional :: mieq
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), optional :: tol
    type (status_t), intent(in out) :: status
    character (*), intent(out) :: msg

    status = NF_STATUS_OK

    if ((present(f_eqcons) .and. .not. present(meq)) .or. &
        .not. present(f_eqcons) .and. present(meq)) then
        msg = "Arguments 'f_eqcons' and 'meq' must be both present"
        goto 100
    end if

    if ((present(f_ieqcons) .and. .not. present(mieq)) .or. &
        .not. present(f_ieqcons) .and. present(mieq)) then
        msg = "Arguments 'f_ieqcons' and 'mieq' must be both present"
        goto 100
    end if

    call check_positive (1, meq, 'meq', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    call check_positive (1, mieq, 'mieq', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    call check_positive (1, maxiter, 'maxiter', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    call check_positive (1.0_PREC, tol, 'tol', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



subroutine slsqp_sanitize_bounds_real64 (xlb_in, xub_in, xlb, xub)
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:), optional :: xlb_in, xub_in
    real (PREC), intent(out), dimension(:) :: xlb, xub

    integer :: n, i
    real (PREC) :: NAN, POS_INF, NEG_INF

    n = size(xlb)

    NAN = ieee_value (0.0_PREC, IEEE_QUIET_NAN)
    POS_INF = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    NEG_INF = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)

    ! Set to NaN by default, as this is the value for signaling that no
    ! lower/upper bound is present in the Scipy-modified code.
    xlb = NAN
    xub = NAN

    if (present(xlb_in)) then
        do i = 1, min(n, size(xlb_in))
            if (ieee_is_finite (xlb_in(i))) then
                xlb(i) = xlb_in(i)
            end if
        end do
    end if

    if (present(xub_in)) then
        do i = 1, min(n, size(xub_in))
            if (ieee_is_finite (xub_in(i))) then
                xub(i) = xub_in(i)
            end if
        end do
    end if

end subroutine



subroutine slsqp_clip_init_real64 (x, xlb, xub)
    integer, parameter :: PREC = real64

    real (PREC), intent(in out), dimension(:) :: x
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
