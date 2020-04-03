

module numfort_optimize_interfaces_real32

    use, intrinsic :: iso_fortran_env

    use numfort_optimize_interfaces_common
    use numfort_common_cond_alloc

    implicit none

    integer, private, parameter :: PREC = real32

    type, public, extends(args_data) :: args_default
        real (PREC), dimension(:), allocatable :: rdata
        integer, dimension(:), allocatable :: idata
    end type

    interface dynamic_cast
        procedure cast_to_args_default
    end interface

    interface cond_alloc
        procedure cond_alloc_args_default
    end interface cond_alloc

    abstract interface
#include "interfaces_spec.F90"
    end interface

    contains

#include "interfaces_impl.F90"

end module
