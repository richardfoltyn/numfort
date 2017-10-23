
#include <numfort.h>

module numfort_optimize_fwrapper

    use, intrinsic :: iso_fortran_env
    use numfort_optimize_interfaces

    implicit none

    private
    public :: dispatch, is_present, wrap_procedure
    public :: fwrapper_ss_real32, fwrapper_ss_real64
    public :: fwrapper_vs_jac_real32, fwrapper_vs_jac_real64
    public :: fwrapper_vv_jac_real32, fwrapper_vv_jac_real64

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
