
module lbfgsb

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: real64
    use numfort_common, only: workspace
    use numfort_optim_result_mod, only: optim_result
    use numfort_optimize_common

    implicit none
    private

    interface minimize_lbfgsb
        module procedure lbfgsb_real64
    end interface

    interface
        function fobj_real64 (x) result(fx)
            import real64
            real (real64), intent(in), dimension(:) :: x
            real (real64) :: fx
        end function

        subroutine grad_real64 (x, fpx)
            import real64
            real (real64), intent(in), dimension(:) :: x
            real (real64), intent(out), dimension(:,:) :: fpx
        end subroutine
    end interface

    public :: workspace, lbfgsb_real64
    public ::
contains

subroutine lbfgsb_real64 (func, x0, grad, lbounds, ubounds, maxiter, maxfun, &
    m, factr, pgtol, iprint, work, res)

    integer, parameter :: PREC = real64

    procedure (fobj_real64), intent(in) :: func
    real (PREC), intent(in out), dimension(:) :: x0
    procedure (grad_real64), intent(in) :: grad
    real (PREC), intent(in), optional :: lbounds, ubounds
    integer, intent(in), optional :: maxiter, maxfun, m, iprint
    real (PREC), intent(in), optional :: factr, pgtol
    class (workspace), intent(in out), optional :: work
    class (optim_result), intent(in out), optional :: res

    real (PREC), dimension(size(x0)) :: llbounds, lubounds
    integer :: lmaxiter, lmaxfun, lm, liprint
    real (PREC) :: lfactr, lpgtol
    type (workspace), pointer :: lwork
    character (len=100) :: msg

    integer :: nrwrk, niwrk, n
    ! lenghts of additional working arrays
    integer, parameter :: ndsave = 29, nisave = 44, nlsave = 4, ncsave = 60

    ! set defaults for optional parameters
    if (present(lbounds)) then
        llbounds = lbounds
    else
        llbounds = IEEE_NEGATIVE_INF
    end if

    if (present(ubounds)) then
        lubounds = ubounds
    else
        lubounds = IEEE_POSITIVE_INF
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

    lfactr = 1d-7
    if (present(factr)) lfactr = factr
    lpgtol = 1d-5
    if (present(pgtol)) lpgtol = pgtol

    ! work array
    n = size(x0)
    niwrk = 3 * n
    nrwrk = 2*lm*n + 5*n + 11*lm*lm + 8*lm

    if (present(work)) then
        lwork => work
    else
        allocate (lwork)
    end if


    return

    ! error handling
100 if (present(res)) then
        res%status = INVALID_INPUT
        res%msg = msg
    end if

    if (.not. present(work)) then
        if (allocated(lwork)) deallocate(lwork)
    end if

end subroutine


end module
