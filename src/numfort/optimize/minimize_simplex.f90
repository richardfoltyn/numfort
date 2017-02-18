! Wrapper for simplex minimizer.
! Author: Richard Foltyn

module numfort_optimize_simplex

    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_optim_result_mod

    use simplex_csiro, only: minim, func_real32, func_real64

    implicit none
    private

    ! add interfaces for 1d wrappers
    interface
        subroutine func_1d_real32(x, fx)
            import real32
            real (real32), intent(in) :: x
            real (real32), intent(out) :: fx
        end subroutine

        subroutine func_1d_real64(x, fx)
            import real64
            real (real64), intent(in) :: x
            real (real64), intent(out) :: fx
        end subroutine
    end interface

    interface minimize_simplex
        module procedure minimize_simplex_real64, minimize_simplex_1d_real64
    end interface

    public :: minimize_simplex

contains

subroutine minimize_simplex_real64 (func, x, tol, maxfun, quad, iprint, res)

    integer, parameter :: PREC = real64
    procedure (func_real64) :: func

    real (PREC), intent(in out), dimension(:) :: x
    real (PREC) :: tol
    integer :: maxfun
    integer (NF_ENUM_KIND) :: iprint
    logical :: quad
    class (optim_result), intent(in out), optional :: res

    intent(in) :: tol, maxfun, quad, iprint
    optional :: maxfun, quad, iprint

    integer :: n, lmaxfun, nloop, iquad, liprint, ifault
    logical :: lquad
    integer (NF_ENUM_KIND) :: status
    real (PREC) :: delta_nonzero, delta_zero, simp, fopt
    real (PREC), dimension(size(x)) :: step, var
    character (100) :: msg

    n = size(x)
    lmaxfun = 20 * n
    lquad = .true.

    if (present(maxfun)) lmaxfun = maxfun
    if (present(quad)) lquad = quad
    iquad = 0
    if (lquad) iquad = 1

    ! disable diagnostic printing
    liprint = -1
    if (present(iprint)) liprint = map_iprint (iprint)

    ! set up initial step size
    ! take some sensible values as used in scipy implementation
    delta_zero = 0.00025d0
    delta_nonzero = 0.05d0
    where (x == 0.0d0)
        step = delta_zero
    else where
        step = x * delta_nonzero
    end where

    nloop = 2 * n
    simp = 1d-6

    call minim (x, step, n, fopt, lmaxfun, liprint, tol, nloop, iquad, &
        simp, var, func, ifault)

    ! map mimum status into numfort_optimize status
    call map_ifault (ifault, status, msg)

    if (status /= NF_STATUS_OK .and. liprint > NF_PRINT_NONE) then
        write (ERROR_UNIT, *) msg
    end if

    if (present(res)) then
        call res%update (x, fopt, status=status, msg=msg)
    end if

end subroutine

subroutine minimize_simplex_1d_real64 (func, x, tol, maxfun, quad, iprint, res)

    integer, parameter :: PREC = real64
    procedure (func_1d_real64) :: func

    real (PREC), intent(in out) :: x
    real (PREC) :: tol
    integer :: maxfun
    integer (NF_ENUM_KIND) :: iprint
    logical :: quad
    class (optim_result), intent(in out), optional :: res

    intent(in) :: tol, maxfun, quad, iprint
    optional :: maxfun, quad, iprint

    real (PREC) :: lx(1)

    lx(1) = x

    call minimize_simplex (func_wrapper, lx, tol, maxfun, quad, iprint, res)

    x = lx(1)

contains
    ! function wrapper that maps calls to objective function with
    ! 1d-array to actual objective with scalar argument
    subroutine func_wrapper (x, fx)
        real (PREC), intent(in), dimension(:) :: x
        real (PREC), intent(out) :: fx

        call func(x(1), fx)
    end subroutine
end subroutine

! Convert numfort_optimize print values to values expected by wrapped routine
pure function map_iprint (i) result(k)
    integer (NF_ENUM_KIND), intent(in) :: i
    integer :: k

    select case (i)
    case (NF_PRINT_MINIMAL)
        k = 0
    case (NF_PRINT_VERBOSE)
        k = 1
    case (NF_PRINT_ALL)
        k = 1
    case default
        ! do not print anything by default
        k = -1
    end select
end function

pure subroutine map_ifault (ifault, status, msg)
    integer, intent(in) :: ifault
        !!  Status code as returned by simplex routine
    integer (NF_ENUM_KIND), intent(out) :: status
        !!  On exit, contains the corresponding NF status code
    character (len=*), intent(out), optional :: msg

    status = NF_STATUS_UNKNOWN
    if (ifault == 0) then
        status = NF_STATUS_OK
        if (present(msg)) msg = "Simplex: successful termination"
    else if (ifault == 1) then
        status = NF_STATUS_MAX_EVAL
        if (present(msg)) msg = "Simplex: max. number of function evaluations exceeded"
    else if (ifault == 3 .or. ifault == 4) then
        status = NF_STATUS_INVALID_ARG
        if (present(msg)) msg = "Simplex: invalid input argument(s)"
    end if
end subroutine

end module
