

module numfort_optimize_slsqp_ng_real64

    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, isfinite => ieee_is_finite

    use numfort_common_enums
    use numfort_common_status
    use numfort_common_input_checks
    use numfort_common_workspace, workspace => workspace_real64

    use numfort_optimize_result_real64
    use numfort_optimize_interfaces_real64
    use numfort_optimize_fwrapper_real64

    use numfort_optimize_slsqp_ng_common

    use numfort_linalg_lh95_real64

    use blas_interfaces, only: GEMV, NRM2

    implicit none

    integer, parameter :: PREC = real64

    private

    public :: minimize_slsqp_ng

    interface minimize_slsqp_ng
        procedure slsqp
    end interface

    integer, parameter :: MAX_INEXACT_LINESEARCH = 10
    integer, parameter :: MAX_RESET_BFGS = 5

    contains

#include "slsqp_ng_impl.f90"

end module
