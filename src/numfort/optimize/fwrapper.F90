
#include <numfort.h>

module numfort_optimize_fwrapper

    use, intrinsic :: iso_fortran_env
    use numfort_optimize_interfaces
    use numfort_optimize_diff

    implicit none

    private
    public :: dispatch
    public :: is_present
    public :: wrap_procedure
    public :: dispatch_fcn_jac

#include <numfort_real32.h>
#include "fwrapper_spec.F90"

#include <numfort_real64.h>
#include "fwrapper_spec.F90"

    contains

#include <numfort_real32.h>
#include "fwrapper_impl.F90"

#include <numfort_real64.h>
#include "fwrapper_impl.F90"


end module
