
#include <numfort.h>


module numfort_polynomial_polyfit

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: vander
    use numfort_common
    use numfort_common_workspace
    use lapack_interfaces, only: GESV

    implicit none
    private

    public :: polyfit
    public :: polyfit_deriv

#include <numfort_real32.h>
#include "polyfit_spec.F90"

#include <numfort_real64.h>
#include "polyfit_spec.F90"

    contains


#include <numfort_real32.h>
#include "polyfit_impl.F90"

#include <numfort_real64.h>
#include "polyfit_impl.F90"

end module
