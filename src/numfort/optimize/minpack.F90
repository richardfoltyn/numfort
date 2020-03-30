module numfort_optimize_minpack

    use, intrinsic :: iso_fortran_env, only: real64

    use numfort_common
    use numfort_common_workspace
    use numfort_common_input_checks
    use numfort_optimize_result
    use numfort_optimize_interfaces
    use numfort_optimize_fwrapper

    use minpack_real64, only: minpack_hybrd_real64 => hybrd, &
        minpack_lmdif_real64 => lmdif, minpack_lmder_real64 => lmder, &
        minpack_hybrj_real64 => hybrj, minpack_chkder => chkder

    implicit none
    private

    public :: root_hybrd
    public :: root_hybrj
    public :: root_lm
    public :: chkder

    integer, parameter :: PREC = real64
    ! size of character variable for diagnostic messages
    integer, parameter :: MSG_LENGTH = 100

    abstract interface
        subroutine func_vec_vec_real64 (x, fx)
            import PREC
            real (PREC), dimension(:), intent(in), contiguous :: x
            real (PREC), dimension(:), intent(out), contiguous :: fx
        end subroutine

        subroutine func_jac_real64 (x, fx, jac, task)
            import PREC
            real (PREC), dimension(:), intent(in), contiguous :: x
            real (PREC), dimension(:), intent(out), contiguous :: fx
            real (PREC), dimension(:,:), intent(out), contiguous :: jac
            integer, intent(inout) :: task
        end subroutine
    end interface

    interface root_hybrd
        module procedure root_hybrd_real64
    end interface

    interface root_hybrj
        module procedure root_hybrj_real64
    end interface

    interface root_lm
        procedure root_lm_real64, &
            root_lm_jac_real64, &
            root_lm_fcn_jac_real64, &
            root_lm_args_real64, &
            root_lm_jac_args_real64, &
            root_lm_fcn_jac_args_real64
    end interface

    interface root_lm_impl
        procedure root_lm_impl_real64
    end interface
    
    interface lm_check_input
        procedure lm_check_input_real64
    end interface

    interface chkder
        module procedure chkder_real64
    end interface

    contains
    
    
#include "minpack_lm_impl.F90"


subroutine root_hybrd_real64 (func, x, fx, xtol, maxfev, ml, mu, eps, factor, diag, work, res)
    integer, parameter :: PREC = real64
    procedure (func_vec_vec_real64) :: func
    ! Note: will be passed using F77 implicit interface, ensure contiguous array.
    real (PREC), dimension(:), contiguous :: x, fx
    real (PREC), dimension(:), target, contiguous :: diag
    real (PREC) :: xtol, eps, factor
    integer :: maxfev, ml, mu
    type (workspace_real64), intent(inout), target, optional :: work
    type (optim_result_real64), intent(inout), optional :: res

    intent (in) :: xtol, maxfev, eps, factor, diag, ml, mu
    intent (inout) :: x, fx
    optional :: xtol, maxfev, eps, factor, diag, ml, mu

    ! local default values for optional arguments
    real (PREC) :: lxtol, leps, lfactor
    integer :: lmaxfev, lml, lmu, lr, lnprint

    ! Remaining local variables
    type (status_t) :: status
    character (MSG_LENGTH) :: msg
    integer :: n, nrwrk, mode, info, i, nfev
    type (workspace_real64), pointer :: ptr_work
    ! pointers to various arrays that need to be passed to hybrd() that are
    ! segments of memory allocated in workspace
    real (PREC), dimension(:), pointer, contiguous :: ptr_fjac, ptr_r, ptr_qtf, &
        ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4, ptr_diag

    nullify (ptr_work)

    n = size(x)
    ! workspace size obtained from hybrd1()
    nrwrk = (n * (3*n + 13)) / 2

    call assert_alloc_ptr (work, ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk)

    ! set default values
    lmaxfev = 200 * (n + 1)
    lml = n - 1
    lmu = n - 1
    lr = n * (n + 1) / 2
    leps = 0.0_PREC
    lfactor = 100.0_PREC
    ! use default from scipy
    lxtol = 1.49012d-08
    lnprint = 0

    mode = 1
    nfev = 0
    info = 0

    ! override defaults if arguments specified by user
    if (present(xtol)) lxtol = xtol
    if (present(ml)) lml = ml
    if (present(mu)) lmu = mu
    if (present(factor)) lfactor = factor
    if (present(eps)) leps = eps

    if (present(diag)) then
        ptr_diag => diag
        ! set mode such that use-provided diag is used
        mode = 2
    else
        ! no scaling of diagonal elements; initialize diag to 1.0 even though
        ! this is not needed
        ptr_diag => ptr_work%rwrk(1:n)
        ptr_diag = 1.0_PREC
    end if

    if (present(maxfev)) then
        if (maxfev > 0) lmaxfev = maxfev
    end if

    ! map array arguments to hybrd into workspace
    ! First n elements reserved for diag
    i = n
    ptr_fjac => ptr_work%rwrk(i+1:i+n*n)
    i = i + n*n
    ptr_r => ptr_work%rwrk(i+1:i+lr)
    i = i + lr
    ptr_qtf => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa1 => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa2 => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa3 => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa4 => ptr_work%rwrk(i+1:i+n)

    call minpack_hybrd_real64 (fwrapper, n, x, fx, lxtol, lmaxfev, lml, lmu, leps, ptr_diag, &
        mode, lfactor, lnprint, info, nfev, ptr_fjac, n, ptr_r, lr, ptr_qtf, &
        ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4)

    select case (info)
    case (0)
        status = NF_STATUS_INVALID_ARG
        msg = "Invalid input parameters"
    case (1)
        status = NF_STATUS_OK
        msg = "Convergence achieved, relative error is at most xtol"
    case (2)
        status = NF_STATUS_MAX_EVAL
        msg = "Exceeded max. number of function evaluations"
    case (3)
        status = NF_STATUS_NOT_CONVERGED
        msg = "xtol too small, no further improvement possible"
    case (4)
        status = NF_STATUS_NOT_CONVERGED
        msg = "Not making progress (as measured by last five Jacobians)"
    case (5)
        status = NF_STATUS_NOT_CONVERGED
        msg = "Not making progress (as measured by last ten iterations)"
    case default
        status = NF_STATUS_UNKNOWN
    end select

    if (present(res)) then
        call result_update (res, x=x, fx=fx, nfev=nfev, nit=1, status=status, msg=msg)
    end if

    call assert_dealloc_ptr (work, ptr_work)

contains

    ! wrapper function: MINPACK's hybrd requires a somewhat inconvenient
    ! function signature, so provide a wrapper around user-supplied function
    subroutine fwrapper (n, x, fx, iflag)
        integer :: n, iflag
        real (PREC) :: x(n), fx(n)

        ! call user-provided function
        call func (x, fx)
    end subroutine

end subroutine


subroutine root_hybrj_real64 (func, x, fx, xtol, maxfev, factor, diag, work, res)
    integer, parameter :: PREC = real64
    procedure (func_jac_real64) :: func
    ! Note: will be passed using F77 implicit interface, ensure contiguous array.
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(inout), dimension(:), contiguous :: fx
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxfev
    real (PREC), intent(in), optional :: factor
    real (PREC), intent(in), dimension(:), target, optional, contiguous :: diag
    type (workspace_real64), intent(inout), target, optional :: work
    type (optim_result_real64), intent(inout), optional :: res

    ! local default values for optional arguments
    real (PREC) :: lxtol, lfactor
    integer :: lmaxfev, lr, lnprint
    ! Remaining local variables
    type (status_t) :: status
    character (MSG_LENGTH) :: msg
    integer :: n, nrwrk, mode, info, i, nfev, njev
    type (workspace_real64), pointer :: ptr_work
    ! pointers to various arrays that need to be passed to hybrd() that are
    ! segments of memory allocated in workspace
    real (PREC), dimension(:), pointer, contiguous :: ptr_fjac, ptr_r, ptr_qtf, &
        ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4, ptr_diag

    nullify (ptr_work)

    n = size(x)
    ! workspace size obtained from hybrj1()
    nrwrk = (n * (3*n + 13)) / 2

    call assert_alloc_ptr (work, ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk)

    ! set default values, taken from hybrj1
    lmaxfev = 100 * (n + 1)
    lr = n * (n + 1) / 2
    lfactor = 100.0_PREC
    ! use default from scipy
    lxtol = 1.49012d-08
    lnprint = 0

    mode = 1
    nfev = 0
    info = 0
    njev = 0

    ! override defaults if arguments specified by user
    if (present(xtol)) lxtol = xtol
    if (present(factor)) lfactor = factor

    if (present(diag)) then
        ptr_diag => diag
        ! set mode such that use-provided diag is used
        mode = 2
    else
        ! no scaling of diagonal elements; initialize diag to 1.0 even though
        ! this is not needed
        ptr_diag => ptr_work%rwrk(1:n)
        ptr_diag = 1.0_PREC
    end if

    if (present(maxfev)) then
        if (maxfev > 0) lmaxfev = maxfev
    end if

    ! map array arguments to hybrj into workspace
    ! First n elements reserved for diag
    i = n
    ptr_fjac => ptr_work%rwrk(i+1:i+n*n)
    i = i + n*n
    ptr_r => ptr_work%rwrk(i+1:i+lr)
    i = i + lr
    ptr_qtf => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa1 => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa2 => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa3 => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa4 => ptr_work%rwrk(i+1:i+n)

    call minpack_hybrj_real64 (fwrapper, n, x, fx, ptr_fjac, n, lxtol, lmaxfev, ptr_diag, &
        mode, lfactor, lnprint, info, nfev, njev, ptr_r, lr, ptr_qtf, &
        ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4)

    select case (info)
    case (0)
        status = NF_STATUS_INVALID_ARG
        msg = "Invalid input parameters"
    case (1)
        status = NF_STATUS_OK
        msg = "Convergence achieved, relative error is at most xtol"
    case (2)
        status = NF_STATUS_MAX_EVAL
        msg = "Exceeded max. number of function evaluations"
    case (3)
        status = NF_STATUS_NOT_CONVERGED
        msg = "xtol too small, no further improvement possible"
    case (4)
        status = NF_STATUS_NOT_CONVERGED
        msg = "Not making progress (as measured by last five Jacobians)"
    case (5)
        status = NF_STATUS_NOT_CONVERGED
        msg = "Not making progress (as measured by last ten iterations)"
    case default
        status = NF_STATUS_UNKNOWN
    end select

    if (present(res)) then
        call result_update (res, x=x, fx=fx, nfev=nfev, nit=1, status=status, msg=msg)
    end if

    call assert_dealloc_ptr (work, ptr_work)

contains

    ! wrapper function: MINPACK's hybrj requires a somewhat inconvenient
    ! function signature, so provide a wrapper around user-supplied function
    subroutine fwrapper (n, x, fx, jac, ldfjac, iflag)
        integer :: n, iflag, ldfjac
        real (PREC) :: x(n), fx(n), jac(ldfjac, n)

        ! call user-provided function
        call func (x, fx, jac, iflag)
    end subroutine

end subroutine



recursive subroutine chkder_real64 (fcn, x, err, status)
    !*  CHKDER verifies that the Jacobian of a function f: R^n -> R^m
    !   is reasonably close to a Jacobian obtained by numerical differentiation.
    procedure (func_jac_real64) :: fcn

    real (PREC), intent(in), dimension(:), contiguous :: x
        !!  Point in function domain where derivative should be evaluated
    real (PREC), intent(out), dimension(:), contiguous :: err
        !!  Array of size m. On exit, indicates which of the 1...m gradients
        !!  of f_i are correct (err(i) = 1.0) or incorrect (err(i) = 0.0).
    type (status_t), intent(out), optional :: status
        !!  Optional status flag.

    real (PREC) :: fvec(size(err)), fvecp(size(err)), fjac(size(err),size(x)), xp(size(x))
    integer :: n, m, mode, task
    type (status_t) :: lstatus

    n = size(x)
    m = size(err)

    mode = 1

    call minpack_chkder (m, n, x, fvec, fjac, m, xp, fvecp, mode, err)

    ! compute fvec, fvecp and jac
    task = 1
    call fcn (xp, fvecp, fjac, task)
    call fcn (x, fvec, fjac, task)
    task = 2
    call fcn (x, fvec, fjac, task)

    mode = 2
    call minpack_chkder (m, n, x, fvec, fjac, m, xp, fvecp, mode, err)

    if (present(status)) status = lstatus

end subroutine

end module
