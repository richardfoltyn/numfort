

module numfort_optimize_simplex_real64
    !*  Wrapper module for CSIRO simplex implementation

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_optimize_result_real64
    use numfort_optimize_interfaces_real64
    use numfort_optimize_fwrapper_real64

    use simplex_csiro, only: minim

    implicit none
    private

    integer, parameter :: PREC = real64

    public :: minimize_simplex

    interface minimize_simplex
        procedure minimize_simplex, minimize_simplex_args
    end interface


    contains

#include "minimize_simplex_impl.F90"

end module
