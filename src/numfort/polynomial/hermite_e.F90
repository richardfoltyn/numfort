

#include <numfort.h>

module numfort_polynomial_hermite_e
    !*  Module contains routines to work with "probabilist" Hermite polynomials,
    !   often abbreviated as H_e.

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_core
    use lapack_interfaces, only: LAPACK_SYEVD => SYEVD

    implicit none
    private

    public :: hermegauss

#include <numfort_real32.h>
#include "herme_spec.F90"

#include <numfort_real64.h>
#include "herme_spec.F90"


    contains


#include <numfort_real32.h>
#include "herme_impl.F90"

#include <numfort_real64.h>
#include "herme_impl.F90"


end module
