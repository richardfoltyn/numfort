
#include "numfort.h"

module numfort_optimize_fminbound
    !*  Module implements derivative-free bounded minimization of scalar
    !   functions.

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_optimize_result
    use numfort_optimize_interfaces
    use numfort_optimize_fwrapper

    implicit none
    private

    public :: minimize_bounded

    interface minimize_bounded
        module procedure minimize_bounded_real32, minimize_bounded_real64
    end interface

    interface minimize_bounded
        module procedure minimize_bounded_args_real32, minimize_bounded_args_real64
    end interface

    interface minimize_bounded_impl
        module procedure minimize_bounded_impl_real32, minimize_bounded_impl_real64
    end interface

    interface check_input
        module procedure check_input_real32, check_input_real64
    end interface

    contains

#define __PREC real32
#include "minimize_bounded_impl.F90"
#undef __PREC

#define __PREC real64
#include "minimize_bounded_impl.F90"
#undef __PREC

end module
