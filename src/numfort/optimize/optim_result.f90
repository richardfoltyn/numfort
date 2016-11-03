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
        real (PREC) :: fx = 0.0_PREC
        real (PREC), dimension(:), allocatable :: x
        integer :: nfev = UNINITIALIZED_COUNTER, nit = UNINITIALIZED_COUNTER
        integer :: status = OPTIM_STATUS_UNKNOWN
        logical :: success = .false.
        character (100) :: msg
    contains
        procedure, pass, private :: update_real64
        procedure, pass, private :: update_real32
        generic, public :: update => update_real32, update_real64

        procedure, pass, private :: reset
    end type

contains

pure subroutine update_real64 (self, x, fx, status, nit, nfev, msg)
    class (optim_result), intent(in out) :: self
    real (real64) :: x(:), fx
    integer :: nit, nfev, status
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg, status
    optional :: x, fx, nit, nfev, msg, status

    integer :: n

    call self%reset ()

    if (present(x)) then
        n = size(x)
        ! allocate array to store optimal point, if needed
        if (allocated(self%x)) then
            if (size(self%x) /= n) then
                deallocate (self%x)
            end if
        end if

        if (.not. allocated(self%x)) then
            allocate (self%x(n))
        end if

        self%x(1:n) = x
    end if

    if (present(fx)) self%fx = fx

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

    call update_real64 (self, x64, real(fx, real64), status, nit, nfev, msg)

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

end module
