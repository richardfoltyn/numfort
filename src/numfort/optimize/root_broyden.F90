
#include <numfort.h>

module numfort_optimize_root_broyden

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_common
    use numfort_common_input_checks
    use numfort_common_workspace

    use numfort_linalg, only: inv, inv_work_query

    use numfort_optimize_result
    use numfort_optimize_interfaces
    use numfort_optimize_fwrapper

    use blas_interfaces, only: DOT, GEMV, SCAL, AXPY, GER
    use lapack_interfaces, only: GETRF, GETRS

    implicit none
    private

    public :: root_broyden

    integer (NF_ENUM_KIND), parameter :: PRINT_LSEARCH = 1
    integer (NF_ENUM_KIND), parameter :: PRINT_STEP = 2
    integer (NF_ENUM_KIND), parameter :: PRINT_JAC = 4

    integer (NF_ENUM_KIND), parameter :: PRINT_NONE = NF_PRINT_NONE
    integer (NF_ENUM_KIND), parameter :: PRINT_ALL = NF_PRINT_ALL

    ! Public fully qualified parameter names
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYDEN_PRINT_LSEARCH = PRINT_LSEARCH
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYDEN_PRINT_STEP = PRINT_STEP
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYDEN_PRINT_JAC = PRINT_JAC
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYDEN_PRINT_ALL = PRINT_ALL
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYEN_PRINT_NONE = PRINT_NONE

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
