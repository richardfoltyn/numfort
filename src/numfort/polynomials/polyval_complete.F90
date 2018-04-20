
#include <numfort.h>

module numfort_polynomial_polyval_complete
    !*  Module implements routine to evaluate complete polynomials,
    !   or their derivatives. This includes routines to compute
    !   complete polynomial basis functions and their Jacobians.

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_core

    use numfort_polynomial_complete

    implicit none
    private

    public :: polybasis_complete
    public :: polybasis_jac_complete

#include <numfort_real32.h>
#include "polyval_complete_spec.F90"

#include <numfort_real64.h>
#include "polyval_complete_spec.F90"

    contains

#include <numfort_real32.h>
#include "polyval_complete_impl.F90"

#include <numfort_real64.h>
#include "polyval_complete_impl.F90"

end module
