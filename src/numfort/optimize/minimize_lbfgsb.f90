
module numfort_optimize_lbfgsb

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: real64
    use numfort_common
    use numfort_optim_result_mod
    use lbfgsb_bmnz_real64, only: setulb

    implicit none
    private

    integer, parameter :: PREC = real64

    integer (NF_ENUM_KIND), parameter :: BOUNDS_NONE = 0
    integer (NF_ENUM_KIND), parameter :: BOUNDS_LOWER = 1
    integer (NF_ENUM_KIND), parameter :: BOUNDS_BOTH = 2
    integer (NF_ENUM_KIND), parameter :: BOUNDS_UPPER = 3

    interface minimize_lbfgsb
        module procedure lbfgsb_real64, lbfgsb_args_real64, lbfgsb_grad_real64, &
            lbfgsb_grad_args_real64
    end interface

    interface
        subroutine fobj_real64 (x, fx)
            import real64
            real (real64), intent(in), dimension(:), contiguous :: x
            real (real64), intent(out) :: fx
        end subroutine

        subroutine grad_real64 (x, fpx)
            import real64
            real (real64), intent(in), dimension(:), contiguous :: x
            real (real64), intent(out), dimension(:), contiguous :: fpx
        end subroutine

        subroutine fobj_grad_real64 (x, fx, fpx)
            import real64
            real (real64), intent(in), dimension(:), contiguous :: x
            real (real64), intent(out) :: fx
            real (real64), intent(out), dimension(:), contiguous :: fpx
        end subroutine

        subroutine fobj_args_real64 (x, fx, args)
            import real64
            real (real64), intent(in), dimension(:), contiguous :: x
            real (real64), intent(out) :: fx
            real (real64), intent(in), dimension(:) :: args
        end subroutine

        subroutine grad_args_real64 (x, fpx, args)
            import real64
            real (real64), intent(in), dimension(:), contiguous :: x
            real (real64), intent(out), dimension(:), contiguous :: fpx
            real (real64), intent(in), dimension(:) :: args
        end subroutine

        subroutine fobj_grad_args_real64 (x, fx, fpx, args)
            import real64
            real (real64), intent(in), dimension(:), contiguous :: x
            real (real64), intent(out) :: fx
            real (real64), intent(out), dimension(:), contiguous :: fpx
            real (real64), intent(in), dimension(:) :: args
        end subroutine
    end interface

    public :: minimize_lbfgsb
contains

! ------------------------------------------------------------------------------
! Wrappers for various ways to invoke the minimizer

subroutine lbfgsb_real64 (func, x, lbounds, ubounds, maxiter, maxfun, &
    m, factr, pgtol, iprint, work, res)

    procedure (fobj_grad_real64) :: func
    real (PREC), intent(in out), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    class (workspace), intent(in out), optional :: work
    class (optim_result), intent(in out), optional :: res

    call lbfgsb_impl_real64 (x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, res=res, func_grad=func)
end subroutine

subroutine lbfgsb_args_real64 (func, x, lbounds, ubounds, maxiter, maxfun, &
    m, factr, pgtol, iprint, work, args, res)

    procedure (fobj_grad_args_real64) :: func
    real (PREC), intent(in out), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    class (workspace), intent(in out), optional :: work
    real (PREC), intent(in), dimension(:) :: args
    class (optim_result), intent(in out), optional :: res

    call lbfgsb_impl_real64 (x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, args=args, res=res, func_grad_args=func)
end subroutine

subroutine lbfgsb_grad_real64 (func, x, grad, lbounds, ubounds, maxiter, maxfun, &
    m, factr, pgtol, iprint, work, res)

    procedure (fobj_real64) :: func
    real (PREC), intent(in out), dimension(:), contiguous :: x
    procedure (grad_real64) :: grad
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    class (workspace), intent(in out), optional :: work
    class (optim_result), intent(in out), optional :: res

    call lbfgsb_impl_real64 (x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, res=res, func=func, grad=grad)
end subroutine

subroutine lbfgsb_grad_args_real64 (func, x, grad, lbounds, ubounds, maxiter, &
    maxfun, m, factr, pgtol, iprint, work, args, res)

    procedure (fobj_args_real64) :: func
    real (PREC), intent(in out), dimension(:), contiguous :: x
    procedure (grad_args_real64) :: grad
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    class (workspace), intent(in out), optional :: work
    real (PREC), intent(in), dimension(:) :: args
    class (optim_result), intent(in out), optional :: res

    call lbfgsb_impl_real64 (x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, args=args, res=res, func_args=func, &
        grad_args=grad)
end subroutine

! ------------------------------------------------------------------------------
! WRAPPER IMPLEMENTATION

subroutine lbfgsb_impl_real64 (x, lbounds, ubounds, maxiter, maxfun, &
        m, factr, pgtol, iprint, work, args, res, &
        func, func_args, grad, grad_args, func_grad, func_grad_args)

    procedure (fobj_real64), optional :: func
    procedure (grad_real64), optional :: grad

    procedure (fobj_args_real64), optional :: func_args
    procedure (grad_args_real64), optional :: grad_args

    procedure (fobj_grad_real64), optional :: func_grad

    procedure (fobj_grad_args_real64), optional :: func_grad_args

    real (PREC), intent(in out), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), optional :: lbounds
    real (PREC), intent(in), dimension(:), optional :: ubounds
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfun
    integer, intent(in), optional :: m
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: factr
    real (PREC), intent(in), optional :: pgtol
    class (workspace), intent(in out), optional, target :: work
    real (PREC), intent(in), dimension(:), optional :: args
    class (optim_result), intent(in out), optional :: res

    real (PREC), dimension(size(x)) :: llbounds, lubounds
    integer :: lmaxiter, lmaxfun, lm
    integer (NF_ENUM_KIND) :: liprint
    integer :: iprint2
        ! iprint argument passed to BFGS library
    real (PREC) :: lfactr, lpgtol
    class (workspace), pointer :: ptr_work
    character (len=100) :: msg

    integer :: nrwrk, niwrk, n
    type (status_t) :: status
    ! lenghts of additional working arrays
    integer, parameter :: ndsave = 29, nisave = 44, nlsave = 4, ncsave = 60, ntask = 60

    ! Allocate all these arrays on the stack since they are small and their
    ! size is independent of the problem.
    character (ncsave) :: csave
    character (ntask) :: task
    logical, dimension(nlsave) :: lsave
    integer, dimension(nisave) :: isave
    real (PREC), dimension(ndsave) :: dsave

    ! lower/upper bound flags passed to setulb
    integer, dimension(size(x)) :: nbd

    ! variables to store function and gradiant
    real (PREC) :: fx, gradx(size(x))
    ! number of iterations and function evaluations
    integer :: iter, nfeval

    status = NF_STATUS_INVALID_ARG
    nullify (ptr_work)

    ! set defaults for optional parameters
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

    lmaxiter = 30
    if (present(maxiter)) lmaxiter = maxiter
    lmaxfun = lmaxiter * 100
    if (present(maxfun)) lmaxfun = maxfun
    lm = 10
    if (present(m)) then
        if (m < 3) then
            msg = 'Number of corrections in limited memory matrix set < 3'
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
        lm = m
    end if

    liprint = NF_PRINT_NONE
    if (present(iprint)) liprint = iprint

    lfactr = 1d+7
    if (present(factr)) lfactr = factr
    lpgtol = 1d-5
    if (present(pgtol)) lpgtol = pgtol

    ! work array dimensions
    n = size(x)
    niwrk = 3 * n
    nrwrk = 2*lm*n + 5*n + 11*lm*lm + 8*lm

    if (present(work)) then
        ptr_work => work
    else
        allocate (workspace :: ptr_work)
    end if

    ! lump together real work array and dsave, and integer work
    ! array and isave, task and csave
    call ptr_work%assert_allocated (nrwrk, niwrk)

    ! transform boundaries to arguments accepted by routine
    call set_bounds_flags (llbounds, lubounds, nbd)
    ! overwrite local print level with value accepted by setulb
    iprint2 = set_print_level (liprint)

    task = 'START'

    iter = 0
    nfeval = 0
    fx = 0.0

    do while (.true.)
        call setulb (n, lm, x, llbounds, lubounds, nbd, fx, gradx, lfactr, lpgtol, &
            ptr_work%rwrk(1:nrwrk), ptr_work%iwrk(1:niwrk), task, iprint2, &
            csave, lsave, isave, dsave)

        if (task(1:2) == 'FG') then
            if (present(func) .and. present(grad)) then
                call func (x, fx)
                call grad (x, gradx)
            else if (present(func_args) .and. present(grad_args) .and. present(args)) then
                call func_args (x, fx, args)
                call grad_args (x, gradx, args)
            else if (present(func_grad)) then
                call func_grad (x, fx, gradx)
            else if (present(func_grad_args) .and. present(args)) then
                call func_grad_args (x, fx, gradx, args)
            else
                msg = "Invalid objective function/gradient specified"
                status = NF_STATUS_INVALID_ARG
                goto 100
            end if
        else if (task(1:5) == 'NEW_X') then
            ! retrieve solver statistics
            iter = isave(30)
            nfeval = isave(34)

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
        msg = task
        status = NF_STATUS_OK
    else if (task(1:5) == 'ERROR') then
        msg = task
        status = NF_STATUS_INVALID_ARG
    else if (task(1:4) == 'ABNO') then
        msg = task
        status = NF_STATUS_NOT_CONVERGED
    else if (iter > lmaxiter) then
        msg = "Number of iterations exceeded limit"
        status = NF_STATUS_MAX_ITER
    else if (nfeval > lmaxfun) then
        msg = "Number of function evaluations exceeded limit"
        status = NF_STATUS_MAX_EVAL
    else
        status = NF_STATUS_UNKNOWN
    end if

100 if (present(res)) then
        ! Note that all variables should have had some value set at this point.
        call res%update (x, fx, status, iter, nfeval, msg)
    end if

    if (.not. present(work) .and. associated(ptr_work)) deallocate(ptr_work)
    nullify (ptr_work)

end subroutine

pure subroutine set_bounds_flags (lb, ub, nbd)
    real (real64), intent(in), dimension(:) :: lb, ub
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

! SET_PRINT_LEVEL converts a verbosity level shared across optimizing
! routines into a value for iprint specific to setulb L-BFGS-B minimizer.
pure function set_print_level (i) result(k)
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
