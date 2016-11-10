module numfort_optimize_minpack

    use, intrinsic :: iso_fortran_env, only: real64

    use numfort_common, only: workspace
    use numfort_optim_result_mod
    use numfort_optimize_common

    use minpack_interfaces, only: hybrd_if

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

    interface root_hybrd
        module procedure root_hybrd_real64
    end interface

    procedure (hybrd_if) :: hybrd

    public :: root_hybrd

contains

subroutine root_hybrd_real64 (f, x, xtol, maxfev, ml, mu, eps, factor, diag, work, res)
    integer, parameter :: PREC = real64
    procedure (func_vec_vec_real64) :: f
    ! Note: will be passed using F77 implicit interface, ensure contiguous array.
    real (PREC), dimension(:), contiguous :: x
    real (PREC), dimension(:), target :: diag
    real (PREC) :: xtol, eps, factor
    integer :: maxfev, ml, mu
    class (workspace), intent(in out), target, optional :: work
    class (optim_result), intent(in out), optional :: res

    intent (in) :: xtol, maxfev, eps, factor, diag, ml, mu
    intent (in out) :: x
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
        ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4, ptr_fvec, ptr_diag

    n = size(x)
    ! on top of what is included in hybrd1, we also place fvec on the workspace
    ! so we need to add n
    nrwrk = (n * (3*n + 13)) / 2 + n

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
        ptr_diag => ptr_work%rwrk(n+1:n+n)
        ptr_diag = 1.0_PREC
    end if

    if (present(maxfev)) then
        if (maxfev > 0) lmaxfev = maxfev
    end if

    ! initial offset: first n elements are reserved for fvec, second segment
    ! if n elements for diag
    ptr_fvec => ptr_work%rwrk(1:n)
    i = 2 * n
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

    call hybrd (fwrapper, n, x, ptr_fvec, lxtol, lmaxfev, lml, lmu, leps, ptr_diag, &
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
        call res%update (x=x, fx=ptr_fvec, nfev=nfev, status=info, msg=ptr_work%cwrk)
    end if

contains

    ! wrapper function: MINPACK's hybrd requires a somewhat inconvenient
    ! function signature, so provide a wrapper around user-supplied function
    subroutine fwrapper (n, x, fx, iflag)
        integer :: n, iflag
        real (PREC) :: x(n), fx(n)

        intent (in) :: n, x
        intent (out) :: fx, iflag

        ! call user-provided function
        call f (x, fx)
    end subroutine

end subroutine

end module
