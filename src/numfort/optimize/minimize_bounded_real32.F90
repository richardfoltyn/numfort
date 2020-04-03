


module numfort_optimize_fminbound_real32
    !*  Module implements derivative-free bounded minimization of scalar
    !   functions.

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_optimize_result_real32
    use numfort_optimize_interfaces_real32
    use numfort_optimize_fwrapper_real32

    implicit none

    integer, parameter :: PREC = real32

    private

    public :: minimize_bounded

    interface minimize_bounded
        procedure minimize_bounded, minimize_bounded_args
    end interface

    contains

#include "minimize_bounded_impl.F90"

end module
