module numfort_optimize_minpack

    use, intrinsic :: iso_fortran_env, only: real64

    use numfort_common, only: workspace
    use numfort_optim_result_mod
    use numfort_optimize_common

    use minpack_interfaces, only: hybrd_if => hybrd, lmdif_if => lmdif

    implicit none
    private

    integer, parameter :: PREC = real64

    interface
        subroutine func_vec_vec_real64 (x, fx)
            import PREC
            real (PREC), dimension(:), intent(in) :: x
            real (PREC), dimension(:), intent(out) :: fx
        end subroutine
    end interface

    interface root_hybr
        module procedure root_hybrd_real64
    end interface

    interface root_lstsq
        module procedure root_lmdif_real64
    end interface

    ! define interfaces for MINPACK procedures
    procedure (hybrd_if) :: hybrd
    procedure (lmdif_if) :: lmdif

    public :: root_hybr, root_lstsq

contains

subroutine root_hybrd_real64 (func, x, fx, xtol, maxfev, ml, mu, eps, factor, diag, work, res)
    integer, parameter :: PREC = real64
    procedure (func_vec_vec_real64) :: func
    ! Note: will be passed using F77 implicit interface, ensure contiguous array.
    real (PREC), dimension(:), contiguous :: x, fx
    real (PREC), dimension(:), target :: diag
    real (PREC) :: xtol, eps, factor
    integer :: maxfev, ml, mu
    class (workspace), intent(in out), target, optional :: work
    class (optim_result), intent(in out), optional :: res

    intent (in) :: xtol, maxfev, eps, factor, diag, ml, mu
    intent (in out) :: x, fx
    optional :: xtol, maxfev, eps, factor, diag, ml, mu

    ! local default values for optional arguments
    real (PREC) :: lxtol, leps, lfactor
    integer :: lmaxfev, lml, lmu, lr, lnprint
    ! Remaining local variables
    integer, parameter :: MSG_LENGTH = 100
    integer :: n, nrwrk, mode, info, i, nfev
    type (workspace), target :: lwork
    class (workspace), pointer :: ptr_work
    ! pointers to various arrays that need to be passed to hybrd() that are
    ! segments of memory allocated in workspace
    real (PREC), dimension(:), pointer, contiguous :: ptr_fjac, ptr_r, ptr_qtf, &
        ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4, ptr_diag

    n = size(x)
    ! workspace size obtained from hybrd1()
    nrwrk = (n * (3*n + 13)) / 2

    if (present(work)) then
        ptr_work => work
    else
        ptr_work => lwork
    end if

    call ptr_work%assert_allocated (nrwrk, ncwrk=MSG_LENGTH)

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

    call hybrd (fwrapper, n, x, fx, lxtol, lmaxfev, lml, lmu, leps, ptr_diag, &
        mode, lfactor, lnprint, info, nfev, ptr_fjac, n, ptr_r, lr, ptr_qtf, &
        ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4)

    select case (info)
    case (0)
        info = OPTIM_STATUS_INVALID_INPUT
        ptr_work%cwrk = "Invalid input parameters"
    case (1)
        info = OPTIM_STATUS_CONVERGED
        ptr_work%cwrk = "Convergence achieved, relative error is at most xtol"
    case (2)
        info = OPTIM_STATUS_MAXFUN
        ptr_work%cwrk = "Exceeded max. number of function evaluations"
    case (3)
        info = OPTIM_STATUS_NOT_CONVERGED
        ptr_work%cwrk = "xtol too small, no further improvement possible"
    case default
        ! convers info = 4 or info = 5
        info = OPTIM_STATUS_NOT_CONVERGED
        ptr_work%cwrk = "Iteration not making good progress"
    end select

    if (present(res)) then
        call res%update (x=x, fx=fx, nfev=nfev, status=info, msg=ptr_work%cwrk)
    end if

contains

    ! wrapper function: MINPACK's hybrd requires a somewhat inconvenient
    ! function signature, so provide a wrapper around user-supplied function
    subroutine fwrapper (n, x, fx, iflag)
        integer :: n, iflag
        real (PREC) :: x(n), fx(n)

        intent (in) :: n, x
        intent (in out) :: fx, iflag

        ! call user-provided function
        call func (x, fx)
    end subroutine

end subroutine

subroutine root_lmdif_real64 (func, x, fx, ftol, xtol, gtol, maxfev, eps, &
    factor, diag, work, res)

    procedure (func_vec_vec_real64) :: func
    ! Note: will be passed using F77 implicit interface, ensure contiguous array.
    real (PREC), dimension(:), contiguous :: x, fx
    real (PREC), dimension(:), target :: diag
    real (PREC) :: xtol, ftol, gtol, eps, factor
    integer :: maxfev
    class (workspace), intent(in out), target, optional :: work
    class (optim_result), intent(in out), optional :: res

    intent (in) :: xtol, ftol, gtol, maxfev, eps, factor, diag
    intent (in out) :: x
    optional :: xtol, ftol, gtol, maxfev, eps, factor, diag

    ! local default values for optional arguments
    real (PREC) :: lxtol, lftol, lgtol, leps, lfactor
    integer :: lmaxfev, lnprint
    ! Remaining local variables
    integer, parameter :: MSG_LENGTH = 100
    integer :: m, n, nrwrk, niwrk, mode, info, i, nfev
    type (workspace), target :: lwork
    class (workspace), pointer :: ptr_work
    ! pointers to various arrays that need to be passed to lmdif() that are
    ! segments of memory allocated in workspace
    real (PREC), dimension(:), pointer, contiguous :: ptr_fjac, ptr_qtf, &
        ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4, ptr_diag
    integer, dimension(:), pointer, contiguous :: ptr_ipvt

    n = size(x)
    m = size(fx)
    ! obtain workspace size from lmdif1()
    nrwrk = m*n + 5*n*m
    ! used to store ipvt
    niwrk = n

    if (present(work)) then
        ptr_work => work
    else
        ptr_work => lwork
    end if
    
    call ptr_work%assert_allocated (nrwrk, niwrk=niwrk, ncwrk=MSG_LENGTH)

    ! set default values
    lmaxfev = 200 * (n + 1)
    leps = 0.0_PREC
    lfactor = 100.0_PREC
    ! use default from scipy
    lxtol = 1.49012d-08
    lftol = 1.49012e-8
    lgtol = 0.0_PREC
    lnprint = 0

    mode = 1
    nfev = 0
    info = 0

    ! override defaults if arguments specified by user
    if (present(xtol)) lxtol = xtol
    if (present(ftol)) lftol = ftol
    if (present(gtol)) lgtol = gtol
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

    ! map array arguments to lmdir() onto workspace
    ! First n elements reserved for diag
    i = n
    ptr_fjac => ptr_work%rwrk(i+1:i+n*m)
    i = i + n*m
    ptr_qtf => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa1 => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa2 => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa3 => ptr_work%rwrk(i+1:i+n)
    i = i + n
    ptr_wa4 => ptr_work%rwrk(i+1:i+m)

    ptr_ipvt => ptr_work%iwrk(1:n)

    call lmdif (fwrapper, m, n, x, fx, lftol, lxtol, lgtol, lmaxfev, &
        leps, ptr_diag, mode, lfactor, lnprint, info, nfev, ptr_fjac, m, &
        ptr_ipvt, ptr_qtf, ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4)

    select case (info)
    case (0)
        info = OPTIM_STATUS_INVALID_INPUT
        ptr_work%cwrk = "Invalid input parameters"
    case (1)
        info = OPTIM_STATUS_CONVERGED
        ptr_work%cwrk = "Convergence in terms of ftol"
    case (2)
        info = OPTIM_STATUS_CONVERGED
        ptr_work%cwrk = "Convergence in terms of xtol"
    case (3)
        info = OPTIM_STATUS_CONVERGED
        ptr_work%cwrk = "Convegence in terms of ftol and xtol"
    case (4)
        info = OPTIM_STATUS_CONVERGED
        ptr_work%cwrk = "Convergence in terms of gtol"
    case (5)
        info = OPTIM_STATUS_MAXFUN
        ptr_work%cwrk = "Exceeded max. number of function evaluations"
    case (6)
        info = OPTIM_STATUS_NOT_CONVERGED
        ptr_work%cwrk = "ftol too small. No further reduction possible"
    case (7)
        info = OPTIM_STATUS_NOT_CONVERGED
        ptr_work%cwrk = "xtol is too small. No further improvement in solution possible"
    case (8)
        ptr_work%cwrk = "gtol is too small. fvec is orthogonal to columns of Jacobian"
    end select

    if (present(res)) then
        call res%update (x=x, fx=fx, nfev=nfev, status=info, msg=ptr_work%cwrk)
    end if


contains

    ! wrapper function: MINPACK's lmdif requires a somewhat inconvenient
    ! function signature, so provide a wrapper around user-supplied function
    subroutine fwrapper (m, n, x, fx, iflag)
        integer :: m, n, iflag
        real (PREC) :: x(n), fx(m)

        intent (in) :: m, n, x
        intent (in out) :: fx, iflag

        ! call user-provided function
        call func (x, fx)
    end subroutine
end subroutine

end module
