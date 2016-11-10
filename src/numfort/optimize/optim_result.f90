! Module implements OPTIM_RESULT, a type used as optimization result
! object for all optimization routines.
! Author: Richard Foltyn

module numfort_optim_result_mod

    use iso_fortran_env, only: real32, real64, int32
    use numfort_common, only: ENUM_KIND

    implicit none
    private

    integer, parameter :: PREC = real64

    integer (ENUM_KIND), parameter :: OPTIM_STATUS_CONVERGED = 0
    integer (ENUM_KIND), parameter :: OPTIM_STATUS_MAXITER = 1
    integer (ENUM_KIND), parameter :: OPTIM_STATUS_MAXFUN = 2
    integer (ENUM_KIND), parameter :: OPTIM_STATUS_NOT_CONVERGED = 4
    integer (ENUM_KIND), parameter :: OPTIM_STATUS_INVALID_INPUT = 2 ** 29
    integer (ENUM_KIND), parameter :: OPTIM_STATUS_UNKNOWN = 2 ** 30

    public :: &
        OPTIM_STATUS_CONVERGED, &
        OPTIM_STATUS_MAXITER, &
        OPTIM_STATUS_MAXFUN, &
        OPTIM_STATUS_NOT_CONVERGED, &
        OPTIM_STATUS_INVALID_INPUT, &
        OPTIM_STATUS_UNKNOWN

    integer, parameter :: UNINITIALIZED_COUNTER = -1

    type, public :: optim_result
        real (PREC), dimension(:), allocatable :: x, fx
        integer :: nfev = UNINITIALIZED_COUNTER, nit = UNINITIALIZED_COUNTER
        integer :: status = OPTIM_STATUS_UNKNOWN
        logical :: success = .false.
        character (100) :: msg
    contains
        procedure, pass, private :: update_real64
        procedure, pass, private :: update_real32
        procedure, pass, private :: update_scalar_scalar_real64
        procedure, pass, private :: update_vec_scalar_real64
        generic, public :: update => update_real32, update_real64, &
            update_scalar_scalar_real64, update_vec_scalar_real64
        procedure, pass, private :: reset
    end type

    interface alloc_assign
        module procedure alloc_assign_real32, alloc_assign_real64
    end interface

contains

! ------------------------------------------------------------------------------
! Implementations for update() generic routine

pure subroutine update_scalar_scalar_real64 (self, x, fx, status, nit, nfev, msg)
    integer, parameter :: PREC = real64
    class (optim_result), intent(in out) :: self
    real (PREC) :: x, fx
    integer :: nit, nfev, status
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg, status
    optional :: nit, nfev, msg, status

    real (PREC), dimension(1) :: x1, fx1

    x1(1) = x
    fx1(1) = fx

    call self%update (x1, fx1, status, nit, nfev, msg)
end subroutine

pure subroutine update_vec_scalar_real64 (self, x, fx, status, nit, nfev, msg)
    integer, parameter :: PREC = real64
    class (optim_result), intent(in out) :: self
    real (PREC) :: x(:), fx
    integer :: nit, nfev, status
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg, status
    optional :: nit, nfev, msg, status

    real (PREC), dimension(1) :: fx1

    fx1(1) = fx

    call self%update (x, fx1, status, nit, nfev, msg)
end subroutine

pure subroutine update_real64 (self, x, fx, status, nit, nfev, msg)
    class (optim_result), intent(in out) :: self
    real (real64), dimension(:) :: x, fx
    integer :: nit, nfev, status
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg, status
    optional :: x, fx, nit, nfev, msg, status

    call self%reset ()

    call alloc_assign (x, self%x)
    call alloc_assign (fx, self%fx)

    if (present(status)) then
        self%success = (status == OPTIM_STATUS_CONVERGED)
        self%status = status
    end if

    if (present(msg)) self%msg = msg

    if (present(nit)) self%nit = nit
    if (present(nfev)) self%nfev = nfev

end subroutine

pure subroutine update_real32 (self, x, fx, status, nit, nfev, msg)
    class (optim_result), intent(in out) :: self
    real (real32) :: x(:), fx
    integer :: nit, nfev, status
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg, status
    ! Note: x cannot be optional in this routine, otherwise generic interface
    ! will not compile
    optional :: fx, nit, nfev, msg, status

    real (real64), dimension(size(x)) :: x64

    x64 = real(x, real64)

    call self%update (x64, real(fx, real64), status, nit, nfev, msg)

end subroutine

pure subroutine reset (self)
    class (optim_result), intent(in out) :: self

    self%nit = UNINITIALIZED_COUNTER
    self%nfev = UNINITIALIZED_COUNTER
    self%msg = ""
    self%fx = 0.0_PREC
    self%status = OPTIM_STATUS_UNKNOWN
    self%success = .false.
end subroutine

pure subroutine alloc_assign_real64 (src, dst)
    integer, parameter :: PREC = real64
    include "include/alloc_assign_impl.f90"
end subroutine

pure subroutine alloc_assign_real32 (src, dst)
    integer, parameter :: PREC = real32
    include "include/alloc_assign_impl.f90"
end subroutine

end module
