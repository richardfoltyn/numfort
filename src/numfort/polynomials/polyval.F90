
#include <numfort.h>

module numfort_polynomial_polyval

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_core

    implicit none
    private

    public :: polyval
    public :: polyder

#include <numfort_real32.h>
#include "polyval_spec.F90"
#include "polyder_spec.F90"

#include <numfort_real64.h>
#include "polyval_spec.F90"
#include "polyder_spec.F90"

    contains

#include <numfort_real32.h>
#include "polyval_impl.F90"
#include "polyder_impl.F90"

#include <numfort_real64.h>
#include "polyval_impl.F90"
#include "polyder_impl.F90"

end module
