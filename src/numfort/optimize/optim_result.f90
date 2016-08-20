

module numfort_optim_result_mod

    use iso_fortran_env

    implicit none

    enum, bind(c)
        enumerator :: OPTIM_STATUS_CONVERGED = 0, &
            OPTIM_STATUS_MAXITER = 1, &
            OPTIM_STATUS_MAXFUN = 2, &
            OPTIM_STATUS_NOT_CONVERGED = 4, &
            OPTIM_STATUS_INVALID_INPUT = 8
    end enum

    type optim_result
        real (real64) :: fx_opt
        real (real64), dimension(:), allocatable :: x_opt
        integer :: nfev = 0, nit = 0, status
        logical :: success = .false.
        character (len=100) :: msg
    end type

end module
