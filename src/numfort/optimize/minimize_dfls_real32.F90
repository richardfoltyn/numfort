

module numfort_optimize_dfls_real32

    use, intrinsic :: iso_fortran_env

    use numfort_common_enums
    use numfort_common_status
    use numfort_common_input_checks
    use numfort_common_workspace, workspace => workspace_real32

    use numfort_optimize_interfaces_real32
    use numfort_optimize_fwrapper_real32
    use numfort_optimize_result_real32

    use newuoa2_real32

    implicit none

    integer, parameter :: PREC = real32

    private

    public :: minimize_dfls

    interface minimize_dfls
        procedure minimize_dfls, minimize_dfls_args
    end interface

    contains

#include "minimize_dfls_impl.F90"

end module
