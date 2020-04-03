

module numfort_optimize_dfls_real64

    use, intrinsic :: iso_fortran_env

    use numfort_common_enums
    use numfort_common_status
    use numfort_common_input_checks
    use numfort_common_workspace, workspace => workspace_real64

    use numfort_optimize_interfaces_real64
    use numfort_optimize_fwrapper_real64
    use numfort_optimize_result_real64

    use newuoa2_real64

    implicit none

    integer, parameter :: PREC = real64

    private

    public :: minimize_dfls

    interface minimize_dfls
        procedure minimize_dfls, minimize_dfls_args
    end interface

    contains

#include "minimize_dfls_impl.F90"

end module
