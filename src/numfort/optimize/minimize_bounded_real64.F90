


module numfort_optimize_fminbound_real64
    !*  Module implements derivative-free bounded minimization of scalar
    !   functions.

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_optimize_result_real64
    use numfort_optimize_interfaces_real64
    use numfort_optimize_fwrapper_real64

    implicit none

    integer, parameter :: PREC = real64

    private

    public :: minimize_bounded

    interface minimize_bounded
        procedure minimize_bounded, minimize_bounded_args
    end interface

    contains

#include "minimize_bounded_impl.F90"

end module
