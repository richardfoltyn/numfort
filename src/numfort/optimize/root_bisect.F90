
#include <numfort.h>

module numfort_optimize_bisect

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_core, only: signum

    use numfort_optimize_result
    use numfort_optimize_interfaces
    use numfort_optimize_fwrapper


    implicit none
    private

    public :: root_bisect

    integer, parameter :: MSG_LEN = 100

#include <numfort_real32.h>
#include "root_bisect_spec.F90"

#include <numfort_real64.h>
#include "root_bisect_spec.F90"

    contains

#include <numfort_real32.h>
#include "root_bisect_impl.F90"

#include <numfort_real64.h>
#include "root_bisect_impl.F90"


end module
