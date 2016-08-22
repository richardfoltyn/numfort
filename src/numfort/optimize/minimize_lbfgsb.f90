
module numfort_optimize_lbfgsb

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: real64
    use numfort_common, only: workspace
    use numfort_optim_result_mod
    use numfort_optimize_common

    implicit none
    private

    integer, parameter :: PREC = real64

    enum, bind(c)
        enumerator :: BOUNDS_NONE = 0, BOUNDS_LOWER = 1, &
            BOUNDS_BOTH = 2, BOUNDS_UPPER = 3
    end enum

    interface minimize_lbfgsb
        module procedure lbfgsb_real64
    end interface

    interface
        subroutine fobj_real64 (x, fx)
            import real64
            real (real64), intent(in), dimension(:) :: x
            real (real64), intent(out) :: fx
        end subroutine

        subroutine grad_real64 (x, fpx)
            import real64
            real (real64), intent(in), dimension(:) :: x
            real (real64), intent(out), dimension(:) :: fpx
        end subroutine
    end interface

    ! interface for F77 setulb, should be useful in eliminating errors at
    ! compile time
    interface
        subroutine if_setulb (n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
            task, iprint, csave, lsave, isave, dsave)

            character (60) :: task, csave
            logical :: lsave(4)
            integer ::n, m, iprint, nbd(n), iwa(3*n), isave(44)
            double precision :: f, factr, pgtol, x(n), l(n), u(n), g(n), &
                wa(2*m*n + 5*n + 11*m*m + 8*m), dsave(29)
        end subroutine
    end interface

    procedure (if_setulb) :: setulb

    public :: minimize_lbfgsb
contains

subroutine lbfgsb_real64 (func, x0, grad, lbounds, ubounds, maxiter, maxfun, &
    m, factr, pgtol, iprint, work, res)

    procedure (fobj_real64) :: func
    real (PREC), intent(in out), dimension(:) :: x0
    procedure (grad_real64) :: grad
    real (PREC), intent(in), dimension(:), optional :: lbounds, ubounds
    integer, intent(in), optional :: maxiter, maxfun, m, iprint
    real (PREC), intent(in), optional :: factr, pgtol
    class (workspace), intent(in out), optional, target :: work
    class (optim_result), intent(in out), optional :: res

    ! actual point in function domain passed to setulb. Need to
    ! have a local copy, otherwise the compiler creates a copy on each
    ! invocation of setulb
    real (PREC), dimension(size(x0)) :: x
    real (PREC), dimension(size(x0)) :: llbounds, lubounds
    integer :: lmaxiter, lmaxfun, lm, liprint
    real (PREC) :: lfactr, lpgtol
    type (workspace), allocatable, target :: lwork
    type (workspace), pointer :: ptr_work => NULL()
    character (len=100) :: msg

    integer :: nrwrk, niwrk, n, status
    ! lenghts of additional working arrays
    integer, parameter :: ndsave = 29, nisave = 44, nlsave = 4, ncsave = 60

    ! working arrays passed to setulb
    real (PREC), dimension(:), pointer, contiguous :: ptr_wa, ptr_dsave
    integer, dimension(:), pointer, contiguous :: ptr_iwa, ptr_isave
    logical, dimension(nlsave) :: lsave
    character (ncsave) :: csave, task

    ! lower/upper bound flags passed to setulb
    integer, dimension(size(x0)) :: nbd

    ! variables to store function and gradiant
    real (PREC) :: fx, gradx(size(x0))
    ! number of iterations and function evaluations
    integer :: iter, nfeval

    status = OPTIM_STATUS_INVALID_INPUT
    x = x0

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
            goto 100
        end if
        lm = m
    end if

    liprint = OPTIM_PRINT_NONE
    if (present(iprint)) liprint = iprint

    lfactr = 1d+7
    if (present(factr)) lfactr = factr
    lpgtol = 1d-5
    if (present(pgtol)) lpgtol = pgtol

    ! work array dimensions
    n = size(x0)
    niwrk = 3 * n
    nrwrk = 2*lm*n + 5*n + 11*lm*lm + 8*lm

    if (present(work)) then
        ptr_work => work
    else
        allocate (lwork)
        ptr_work => lwork
    end if

    ! lump together real work array and dsave, and integer work
    ! array and isave
    call ptr_work%assert_allocated (nrwrk + ndsave, niwrk + nisave, ncsave, nlsave)

    ptr_wa => ptr_work%rwrk(1:nrwrk)
    ptr_iwa => ptr_work%iwrk(1:niwrk)
    ptr_dsave => ptr_work%rwrk(nrwrk+1:)
    ptr_isave => ptr_work%iwrk(niwrk+1:)

    ! transform boundaries to arguments accepted by routine
    call set_bounds_flags (llbounds, lubounds, nbd)
    ! overwrite local print level with value accepted by setulb
    liprint = set_print_level (liprint)

    task = 'START'

    iter = 0
    nfeval = 0

    do while (.true.)
        call setulb (n, lm, x, llbounds, lubounds, nbd, fx, gradx, lfactr, lpgtol, &
            ptr_wa, ptr_iwa, task, liprint, csave, lsave, ptr_isave, ptr_dsave)

        if (task(1:2) == 'FG') then
            call func (x, fx)
            call grad (x, gradx)
        else if (task(1:5) == 'NEW_X') then
            ! retrieve statitics stored in working array
            iter = ptr_isave(30)
            nfeval = ptr_isave(34)
            if (iter > lmaxiter .or. nfeval > lmaxfun) exit
        else
            exit
        end if
    end do

    if (task(1:4) == 'CONV') then
        msg = task
        status = OPTIM_STATUS_CONVERGED
    else if (task(1:5) == 'ERROR') then
        msg = task
        status = OPTIM_STATUS_INVALID_INPUT
    else if (task(1:4) == 'ABNO') then
        msg = task
        status = OPTIM_STATUS_NOT_CONVERGED
    else if (iter > lmaxiter) then
        msg = "STOP: number of iterations exceeded limit"
        status = OPTIM_STATUS_MAXITER
    else if (nfeval > lmaxfun) then
        msg = "STOP: number of function evaluations exceeded limit"
        status = OPTIM_STATUS_MAXFUN
    end if

    ! store found minimum back into argument array
    x0 = x

    if (present(res)) then
        call res%update (x, fx, status, iter, nfeval, msg)
    end if

    return

    ! input error handling
100 if (present(res)) then
        res%status = OPTIM_STATUS_INVALID_INPUT
        res%msg = msg
    end if

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
    integer, intent(in) :: i
    integer :: k

    select case (i)
    case (OPTIM_PRINT_MINIMAL)
        k = 10
    case (OPTIM_PRINT_VERBOSE)
        k = 99
    case (OPTIM_PRINT_ALL)
        k = 101
    case default
        ! do not print anything by default
        k = -1
    end select
end function

end module
