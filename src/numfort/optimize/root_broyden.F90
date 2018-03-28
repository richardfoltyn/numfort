
#include <numfort.h>

module numfort_optimize_root_broyden

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_common_workspace

    use numfort_linalg, only: inv

    use numfort_optimize_result
    use numfort_optimize_interfaces
    use numfort_optimize_fwrapper

    use blas_interfaces, only: DOT, GEMV, SCAL, AXPY, GER
    use lapack_interfaces, only: GETRF, GETRS


    implicit none
    private

    public :: root_broyden

    integer, parameter :: LINESEARCH_MAX_STEPS = 4
        !*  Max. number of evaluations during line search (ie. at most
        !   LINESEARCH_MAX_STEPS-1 backtracking steps will be performed)

#include <numfort_real32.h>
#include "root_broyden_spec.F90"

#include <numfort_real64.h>
#include "root_broyden_spec.F90"

    contains


#include <numfort_real32.h>
#include "root_broyden_impl.F90"

#include <numfort_real64.h>
#include "root_broyden_impl.F90"

end module
