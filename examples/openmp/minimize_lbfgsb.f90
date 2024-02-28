! Demo of parallelized L-BFGS-B use in a very simplified economic example:
! For a set of asset levels a_1,...,a_m, the following maximization problem
! is solved in a parallelized fashion:
!   V_t(a) = max_a' {u(a-a') + beta * V_{t+1}(a')}
! for only one period t.
! For this purpose, it is assumed that V_{t+1} = u(a') / (1-beta)

module problem_mod

    use, intrinsic :: iso_fortran_env
    use numfort_optimize, workspace => workspace_real64, &
        optim_result => optim_result_real64
    implicit none
    private

    public :: solve_savings

    integer, parameter :: PREC = real64
    real (PREC), parameter :: sigma = 2.0d0

contains

subroutine solve_savings (a_at, beta, work, sav)
    real (PREC), intent(in) :: a_at, beta
    type (workspace), intent(in out) :: work
    real (PREC), intent(out) :: sav

    real (PREC) :: x(1), lbnd(1), ubnd(1), args(2)
    integer (NF_ENUM_KIND) :: iprint
    type (optim_result) :: res

    iprint = NF_PRINT_NONE

    ! Initial guess: save half of asset holdings
    sav = a_at / 2.0d0
    ! Bounds on feasible savings
    lbnd(1) = 0.0_PREC
    ubnd(1) = a_at
    ! store some additional data that we want to pass on to objective function
    args(1) = a_at
    args(2) = beta

    x(1) = sav
    call minimize_lbfgsb (fobj_grad, x, lbounds=lbnd, ubounds=ubnd, &
        iprint=iprint, work=work, res=res)
    sav = x(1)

    contains

    subroutine fobj_grad (x, fx, g)
        real (PREC), intent(in), dimension(:), contiguous :: x
        real (PREC), intent(out) :: fx
        real (PREC), intent(out), dimension(:), contiguous :: g

        real (PREC) :: cons, vcont, ap, a, sav

        sav = x(1)
        a = a_at

        cons = a - sav
        ! Assume that interest rate is given by r = 1/beta - 1
        ap = sav / beta
        ! Assume that continuation value is given by u(a') / (1-beta)
        vcont = util (ap) / (1.0d0 - beta)
        ! Objective is formulated in terms of minimization
        fx = - (util(cons) + beta * vcont)
        g(1) = util_c(cons) - util_c(ap) / (1.0d0 - beta)
    end subroutine

end subroutine

pure function util (c) result(res)
    real (PREC), intent(in) :: c
    real (PREC) :: res

    if (c <= 0.0) then
        res = -1.0d20
    else
        res = (c ** (1.0d0 - sigma) - 1) / (1.0d0 - sigma)
    end if
end function

pure function util_c (c) result(res)
    real (PREC), intent(in) :: c
    real (PREC) :: res

    if (c <= 0.0) then
        res = 1.0d20
    else
        res = c ** (-sigma)
    end if
end function

end module


program lbfgsb_omp

    use numfort_optimize, workspace => workspace_real64
    use iso_fortran_env
    use omp_lib
    use problem_mod

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()

contains


subroutine example1 ()

    real (PREC), parameter :: beta = 0.95d0
    integer, parameter :: m = int(1e4)
    real (PREC), parameter :: amin = 0.1d0, amax = 100.0d0

    ! asset grid
    real (PREC), dimension(:), allocatable :: agrid
    real (PREC), dimension(:), allocatable :: sav_opt

    ! Variables private to threads
    type (workspace) :: ws
    real (PREC) :: s, a_at
    integer :: i

    allocate (agrid(m), sav_opt(m))

    ! create power-spaced asset grid
    call power_grid (agrid, amin, amax, 3)

    !$omp parallel default(none) &
    !$omp private(i, ws, s, a_at) &
    !$omp shared(agrid, sav_opt)

    !$omp do schedule(auto)
    do i = 1, m
        a_at = agrid(i)
        call solve_savings (a_at, beta, ws, s)
        sav_opt(i) = s
    end do
    !$omp end do

    !$omp end parallel

    ! report results
    print *, "Optimal savings:"
    print '(2(tr2, a10, tr3, a10))', "Assets", "Savings", "Assets", "Savings"
    print '(*(2(tr2, f10.6, tr3, f10.6), /))', (agrid(i), sav_opt(i), i=1,m)

end subroutine

! POWER_GRID creates a "power-spaced" grid from xmin to xmax
subroutine power_grid(x, xmin, xmax, pow)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:) :: x
    real (PREC), intent(in) :: xmin, xmax
    integer, intent(in) :: pow

    real (real64), dimension(size(x)) :: work
    integer :: i, n

    n = size(x)

    do i = 1, n 
        work(i) = ((i-1.0d0) / (n - 1.0d0)) ** pow
    end do
    x = real(xmin + (xmax - xmin) * work, PREC)
    ! prevent any rounding errors
    x(1) = xmin
    x(n) = xmax

end subroutine

end program
