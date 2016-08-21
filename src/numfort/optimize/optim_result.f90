

module numfort_optim_result_mod

    use iso_fortran_env, only: real32, real64, int32
    use numfort_common, only: ENUM_KIND

    implicit none
    private

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

    type, public :: optim_result
        real (real64) :: fx_opt
        real (real64), dimension(:), allocatable :: x_opt
        integer :: nfev = -1, nit = -1, status = OPTIM_STATUS_UNKNOWN
        logical :: success = .false.
        character (len=100) :: msg
    contains
        procedure, pass, private :: update_real64
        procedure, pass, private :: update_real32
        generic, public :: update => update_real32, update_real64
    end type

contains

pure subroutine update_real64 (self, x, fx, status, nit, nfev, msg)
    class (optim_result), intent(in out) :: self
    real (real64) :: x(:), fx
    integer :: nit, nfev, status
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg, status
    optional :: fx, nit, nfev, msg, status

    integer :: n

    n = size(x)

    if (present(fx)) self%fx_opt = fx
    if (present(status)) then
        self%success = (status == OPTIM_STATUS_CONVERGED)
        self%status = status
    end if

    if (present(msg)) self%msg = msg

    ! allocate array to store minimum, if needed
    if (allocated(self%x_opt)) then
        if (size(self%x_opt) < n) then
            deallocate (self%x_opt)
            allocate (self%x_opt(n), source=x)
        else
            self%x_opt(1:n) = x
        end if
    else
        allocate (self%x_opt(n), source=x)
    end if

    if (present(nit)) self%nit = nit
    if (present(nfev)) self%nfev = nfev

end subroutine

pure subroutine update_real32 (self, x, fx, status, nit, nfev, msg)
    class (optim_result), intent(in out) :: self
    real (real32) :: x(:), fx
    integer :: nit, nfev, status
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg, status
    optional :: fx, nit, nfev, msg, status

    real (real64), dimension(size(x)) :: x64

    x64 = real(x, real64)

    call update_real64 (self, x64, real(fx, real64), status, nit, nfev, msg)

end subroutine

end module
