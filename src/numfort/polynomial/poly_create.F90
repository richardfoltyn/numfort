

#include <numfort.h>


module numfort_polynomial_create
    !*  Module implements various routines to create polynomials that do
    !   not fit into a more specific module (such as fitting polynomials).

    use, intrinsic :: iso_fortran_env

    use numfort_core
    use numfort_common

    implicit none
    private

    public :: polyder
    public :: polyshift

#include <numfort_real32.h>
#include "poly_create_spec.F90"

#include <numfort_real64.h>
#include "poly_create_spec.F90"

    contains

#include <numfort_real32.h>
#include "poly_create_impl.F90"

#include <numfort_real64.h>
#include "poly_create_impl.F90"


end module