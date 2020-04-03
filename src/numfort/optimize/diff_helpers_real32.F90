
module numfort_optimize_diff_real32
    !*  Module contains helper routines to numerically differentiate
    !   (multivariate) functions.

    use, intrinsic :: iso_fortran_env

    use numfort_optimize_interfaces_real32

    implicit none

    private
    public :: num_diff

    integer, parameter :: PREC = real32

    interface num_diff
        procedure fss_deriv, fss_deriv_args, fvs_deriv, fvs_args_deriv, &
            fvv_deriv, fvv_args_deriv
    end interface

    interface get_step_size
        procedure get_step_size
    end interface


    contains

#include "diff_helpers_impl.F90"

end module
