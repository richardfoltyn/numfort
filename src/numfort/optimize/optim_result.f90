

module numfort_optim_result_mod

    use iso_fortran_env

    implicit none

    enum, bind(c)
        enumerator :: OPTIM_STATUS_INVALID_INPUT = 100
    end enum

    type optim_result
        real (real64) :: fx_opt
        real (real64), dimension(:), allocatable :: x_opt
        integer :: feval, iter, status
        character (len=100) :: msg
    end type

end module
