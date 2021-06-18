

module numfort_arrays_grid_real64
    !*  Module implements grid creation routines using various spacing
    !   algorithms.

    use, intrinsic :: iso_fortran_env

    use numfort_core
    use numfort_common_status
    use numfort_common_input_checks

    implicit none

    private

    public :: linspace
    public :: powerspace
    public :: logspace
    public :: log_shift_space

    integer, parameter :: PREC = real64

    interface linspace
        procedure linspace
    end interface

    interface powerspace
        procedure powerspace
    end interface

    interface logspace
        procedure logspace
    end interface

    interface log_shift_space
        procedure log_shift_space
    end interface

    contains

#include "grid_impl.F90"


end module
