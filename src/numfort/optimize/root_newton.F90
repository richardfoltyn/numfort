
#include <numfort.h>

module numfort_optimize_newton

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_core, only: signum

    use numfort_optimize_result
    use numfort_optimize_interfaces
    use numfort_optimize_fwrapper


    implicit none
    private

    public :: root_newton
    public :: root_halley

    integer, parameter :: MSG_LEN = 100

#include <numfort_real32.h>
#include "root_newton_spec.F90"

#include <numfort_real64.h>
#include "root_newton_spec.F90"

    contains

#include <numfort_real32.h>
#include "root_newton_impl.F90"

#include <numfort_real64.h>
#include "root_newton_impl.F90"


end module
