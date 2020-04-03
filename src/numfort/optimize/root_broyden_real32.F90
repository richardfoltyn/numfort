

module numfort_optimize_root_broyden_real32

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_common
    use numfort_common_input_checks
    use numfort_common_workspace, workspace => workspace_real32

    use numfort_linalg, only: inv, inv_work_query

    use numfort_optimize_result_real32
    use numfort_optimize_interfaces_real32
    use numfort_optimize_fwrapper_real32
    use numfort_optimize_root_broyden_common

    use blas_interfaces, only: COPY, DOT, GEMV, SCAL, AXPY, GER
    use lapack_interfaces, only: GETRF, GETRS

    implicit none

    integer, parameter :: PREC = real32

    private

    public :: root_broyden

    interface root_broyden
        procedure root_broyden, root_broyden_args, &
            root_broyden_jac, root_broyden_jac_args, &
            root_broyden_fcn_jac_opt, root_broyden_fcn_jac_opt_args
    end interface

    contains

#include "root_broyden_impl.F90"

end module
