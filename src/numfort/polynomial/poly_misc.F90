

#include <numfort.h>

module numfort_polynomial_misc
    !*  Module implements misc. polynomial (helper) functions that do
    !   not fit in any more specific polynomial module.

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_core_combinatorics

    public :: polyshift

#include <numfort_real32.h>
#include "poly_misc_spec.F90"

#include <numfort_real64.h>
#include "poly_misc_spec.F90"

    contains

#include <numfort_real32.h>
#include "poly_misc_impl.F90"

#include <numfort_real64.h>
#include "poly_misc_impl.F90"

end module