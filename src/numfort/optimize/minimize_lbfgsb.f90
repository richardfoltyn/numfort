
module numfort_optimize_lbfgsb

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_common_workspace

    use numfort_optimize_result
    use numfort_optimize_fwrapper
    use numfort_optimize_interfaces

    use lbfgsb_bmnz_real64, only: setulb

    implicit none
    private

    public :: minimize_lbfgsb

    integer (NF_ENUM_KIND), parameter :: BOUNDS_NONE = 0
    integer (NF_ENUM_KIND), parameter :: BOUNDS_LOWER = 1
    integer (NF_ENUM_KIND), parameter :: BOUNDS_BOTH = 2
    integer (NF_ENUM_KIND), parameter :: BOUNDS_UPPER = 3

    interface minimize_lbfgsb
        procedure lbfgsb_real64, lbfgsb_args_real64, &
            lbfgsb_jac_real64, lbfgsb_jac_args_real64, &
            lbfgsb_fcn_jac_real64, lbfgsb_fcn_jac_args_real64
    end interface

    interface set_bounds_flags
       procedure set_bounds_flags_real64
    end interface

    interface check_inputs
        procedure check_inputs_real64
    end interface

    contains


subroutine check_inputs_real64 (maxiter, maxfun, m, factr, pgtol, res)
    integer, parameter :: PREC = real64
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    type (optim_result_real64), intent(inout) :: res

    res%status = NF_STATUS_OK

    call check_positive (0.0_PREC, factr, "factr", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, pgtol, "pgtol", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (0, m, "m", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (0, maxiter, "maxiter", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    call check_positive (0, maxfun, "maxfun", res%status, res%msg)
    if (res%status /= NF_STATUS_OK) goto 100

    if (present(m)) then
        if (m < 3) then
            res%msg = 'Number of corrections in limited memory matrix set < 3'
            res%status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    return

100 continue
    res%status = NF_STATUS_INVALID_ARG

end subroutine

! ------------------------------------------------------------------------------
! Wrappers for various ways to invoke the minimizer

recursive subroutine lbfgsb_real64 (fcn, x, ndiff, lbounds, ubounds, maxiter, &
        maxfun, m, factr, pgtol, iprint, dstep, work, res)

    integer, parameter :: PREC = real64
    procedure (fvs_fcn_real64) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    logical, intent(in) :: ndiff
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    real (PREC), intent(in), optional :: dstep
        !*  If present, sets the step size for numerical differentiation.
    type (workspace_real64), intent(inout), optional :: work
    type (optim_result_real64), intent(inout), optional :: res

    type (fwrapper_vs_real64) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            call result_reset (res)
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn=fcn, eps=dstep)

    call lbfgsb_impl_real64 (fwrapper, x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, res)
end subroutine


recursive subroutine lbfgsb_args_real64 (fcn, x, args, ndiff, lbounds, ubounds, &
        maxiter, maxfun, m, factr, pgtol, iprint, dstep, work, res)

    integer, parameter :: PREC = real64
    procedure (fvs_fcn_real64) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    logical, intent(in) :: ndiff
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    real (PREC), intent(in), optional :: dstep
    !*  If present, sets the step size for numerical differentiation.
    type (workspace_real64), intent(inout), optional :: work
    type (optim_result_real64), intent(inout), optional :: res

    type (fwrapper_vs_real64) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            call result_reset (res)
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn=fcn, args=args, eps=dstep)

    call lbfgsb_impl_real64 (fwrapper, x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, res)
end subroutine


recursive subroutine lbfgsb_jac_real64 (fcn, fjac, x, lbounds, ubounds, &
        maxiter, maxfun, m, factr, pgtol, iprint, work, res)

    integer, parameter :: PREC = real64
    procedure (fvs_fcn_real64) :: fcn
    procedure (fvs_jac_real64) :: fjac
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    type (workspace_real64), intent(inout), optional :: work
    type (optim_result_real64), intent(inout), optional :: res

    type (fwrapper_vs_real64) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn, jac=fjac)

    call lbfgsb_impl_real64 (fwrapper, x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, res)
end subroutine


recursive subroutine lbfgsb_jac_args_real64 (fcn, fjac, x, args, &
        lbounds, ubounds, maxiter, maxfun, m, factr, pgtol, iprint, work, res)

    integer, parameter :: PREC = real64
    procedure (fvs_fcn_real64) :: fcn
    procedure (fvs_jac_real64) :: fjac
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    type (workspace_real64), intent(inout), optional :: work
    type (optim_result_real64), intent(inout), optional :: res

    type (fwrapper_vs_real64) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn, jac=fjac, args=args)

    call lbfgsb_impl_real64 (fwrapper, x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, res)
end subroutine


recursive subroutine lbfgsb_fcn_jac_real64 (fcn, x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, res)

    integer, parameter :: PREC = real64
    procedure (fvs_fcn_jac_real64) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    type (workspace_real64), intent(inout), optional :: work
    type (optim_result_real64), intent(inout), optional :: res

    type (fwrapper_vs_real64) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac=fcn)

    call lbfgsb_impl_real64 (fwrapper, x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, res)
end subroutine

recursive subroutine lbfgsb_fcn_jac_args_real64 (fcn, x, args, lbounds, ubounds,&
        maxiter, maxfun, m, factr, pgtol, iprint, work, res)

    integer, parameter :: PREC = real64
    procedure (fvs_fcn_jac_real64) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    type (workspace_real64), intent(inout), optional :: work
    type (optim_result_real64), intent(inout), optional :: res

    type (fwrapper_vs_real64) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac=fcn, args=args)

    call lbfgsb_impl_real64 (fwrapper, x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, res)
end subroutine


! ------------------------------------------------------------------------------
! WRAPPER IMPLEMENTATION

recursive subroutine lbfgsb_impl_real64 (fcn, x, lbounds, ubounds, maxiter, &
        maxfun, m, factr, pgtol, iprint, work, res)

    integer, parameter :: PREC = real64

    type (fwrapper_vs_real64) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    type (workspace_real64), intent(inout), optional, target :: work
    type (optim_result_real64), intent(inout), optional, target :: res

    type (optim_result_real64), pointer :: ptr_res

    integer :: lmaxiter, lmaxfun, lm
    integer (NF_ENUM_KIND) :: liprint
    integer :: iprint2
        ! iprint argument passed to BFGS library
    real (PREC) :: lfactr, lpgtol
    type (workspace_real64), pointer :: ptr_work

    integer :: nrwrk, niwrk, n, lwork, liwork
    type (status_t) :: status
    ! lenghts of additional working arrays
    integer, parameter :: NDSAVE = 29, NISAVE = 44, NLSAVE = 4, NCSAVE = 60, NTASK = 60

    ! Allocate all these arrays on the stack since they are small and their
    ! size is independent of the problem.
    character (NCSAVE) :: csave
    character (NTASK) :: task
    logical, dimension(NLSAVE) :: lsave

    real (PREC) :: fx
    integer, dimension(:), pointer, contiguous :: nbd
        !  lower/upper bound flags passed to setulb
    real (PREC), dimension(:), pointer, contiguous :: fpx
        !   Array to store gradient
    real (PREC), dimension(:), pointer, contiguous :: ptr_rwork
    integer, dimension(:), pointer, contiguous :: ptr_iwork
    real (PREC), dimension(:), pointer, contiguous :: llbounds, lubounds
        !   Local arrays storing lower/upper bounds
    real (PREC), dimension(:), pointer, contiguous :: ptr_dsave
    integer, dimension(:), pointer, contiguous :: ptr_isave

    integer :: iter, nfeval
        !   Number of iterations and function evaluations

    status = NF_STATUS_INVALID_ARG
    nullify (ptr_work)
    nullify (ptr_rwork, ptr_iwork, llbounds, lubounds, fpx)

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call check_inputs (maxiter, maxfun, m, factr, pgtol, ptr_res)
    if (ptr_res%status /= NF_STATUS_OK) goto 100

    ! Default parameter values
    lmaxiter = 30
    lmaxfun = lmaxiter * 100
    lm = 10
    liprint = NF_PRINT_NONE
    lfactr = 1d+7
    lpgtol = 1d-5

    if (present(maxiter)) lmaxiter = maxiter
    if (present(maxfun)) lmaxfun = maxfun
    if (present(iprint)) liprint = iprint
    if (present(factr)) lfactr = factr
    if (present(pgtol)) lpgtol = pgtol

    ! work array dimensions
    n = size(x)
    liwork = 3 * n
    lwork = 2*lm*n + 5*n + 11*lm*lm + 8*lm

    ! Add space for bounds, fpx
    nrwrk = lwork + 3*n
    ! Add space for additional working arrays used by SETULB
    nrwrk = nrwrk + NDSAVE
    ! Add space for NBD
    niwrk = liwork + n
    ! Additional working arrays
    niwrk = niwrk * NISAVE

    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal workspace array offsets
    call workspace_reset (ptr_work)

    call assert_alloc (ptr_work, nrwrk=nrwrk, niwrk=niwrk)

    ! Initialize pointers to workspace arrays
    call workspace_get_ptr (ptr_work, n, fpx)
    call workspace_get_ptr (ptr_work, n, llbounds)
    call workspace_get_ptr (ptr_work, n, lubounds)
    call workspace_get_ptr (ptr_work, n, nbd)

    call workspace_get_ptr (ptr_work, lwork, ptr_rwork)
    call workspace_get_ptr (ptr_work, liwork, ptr_iwork)

    call workspace_get_ptr (ptr_work, NDSAVE, ptr_dsave)
    call workspace_get_ptr (ptr_work, NISAVE, ptr_isave)

    ! Set lower / upper bounds
    if (present(lbounds)) then
        llbounds = lbounds
    else
        llbounds = ieee_value (llbounds(1), IEEE_NEGATIVE_INF)
    end if

    if (present(ubounds)) then
        lubounds = ubounds
    else
        lubounds = ieee_value (lubounds(1), IEEE_POSITIVE_INF)
    end if

    ! transform boundaries to arguments accepted by routine
    call set_bounds_flags (llbounds, lubounds, nbd)
    ! overwrite local print level with value accepted by setulb
    iprint2 = set_print_level (liprint)

    task = 'START'

    iter = 0
    nfeval = 0
    fx = 0.0

    do while (.true.)
        call setulb (n, lm, x, llbounds, lubounds, nbd, fx, fpx, lfactr, lpgtol, &
            ptr_rwork, ptr_iwork, task, iprint2, csave, lsave, ptr_isave, ptr_dsave)

        if (task(1:2) == 'FG') then
            call dispatch_fcn_jac (fcn, x, fx=fx, fpx=fpx)
        else if (task(1:5) == 'NEW_X') then
            ! retrieve solver statistics
            iter = ptr_isave(30)
            nfeval = ptr_isave(34)

            ! check whether limits have been exceeded
            if (iter > lmaxiter) then
                task = 'STOP: Number of iterations exceeds limit'
            else if (nfeval > lmaxfun) then
                task = 'STOP: Total number of function evaluations exceeds limit'
            end if
        else
            exit
        end if
    end do

    if (task(1:4) == 'CONV') then
        ptr_res%msg = task
        ptr_res%status = NF_STATUS_OK
    else if (task(1:5) == 'ERROR') then
        ptr_res%msg = task
        ptr_res%status = NF_STATUS_INVALID_ARG
    else if (task(1:4) == 'ABNO') then
        ptr_res%msg = task
        ptr_res%status = NF_STATUS_NOT_CONVERGED
    else if (iter > lmaxiter) then
        ptr_res%msg = "Number of iterations exceeded limit"
        ptr_res%status = NF_STATUS_MAX_ITER
    else if (nfeval > lmaxfun) then
        ptr_res%msg = "Number of function evaluations exceeded limit"
        ptr_res%status = NF_STATUS_MAX_EVAL
    else
        ptr_res%status = NF_STATUS_UNKNOWN
        ptr_res%msg = "Unknown exit code"
    end if

100 continue

    call result_update (ptr_res, x, fx, nit=iter, nfev=nfeval)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

end subroutine


pure subroutine set_bounds_flags_real64 (lb, ub, nbd)
    !*  SET_BOUNDS_FLAG sets the boundary flag expected by the underlying
    !   L-BFGS-B implementation depending on user-provided bounds
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:) :: lb, ub
    integer, dimension(:), intent(out) :: nbd

    integer :: i

    do i = 1, size(nbd)
        if (ieee_is_finite(lb(i)) .and. ieee_is_finite(ub(i))) then
            nbd(i) = BOUNDS_BOTH
        else if (ieee_is_finite(lb(i))) then
            nbd(i) = BOUNDS_LOWER
        else if (ieee_is_finite(ub(i))) then
            nbd(i) = BOUNDS_UPPER
        else
            nbd(i) = BOUNDS_NONE
        end if
    end do
end subroutine


pure function set_print_level (i) result(k)
    !*  SET_PRINT_LEVEL converts a verbosity level shared across optimizing
    !   routines into a value for iprint specific to setulb L-BFGS-B minimizer.
    integer (NF_ENUM_KIND), intent(in) :: i
    integer :: k

    select case (i)
    case (NF_PRINT_MINIMAL)
        k = 10
    case (NF_PRINT_VERBOSE)
        k = 99
    case (NF_PRINT_ALL)
        k = 101
    case default
        ! do not print anything by default
        k = -1
    end select
end function

end module
