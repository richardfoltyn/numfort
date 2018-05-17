
#include <numfort.h>

module numfort_polynomial_ppolyval

    use, intrinsic :: iso_fortran_env

    use numfort_core
    use numfort_common
    use numfort_common_status
    use numfort_common_workspace
    use numfort_interpolate_common
    use numfort_polynomial_ppoly

    implicit none
    private

    public :: ppolyval

#include <numfort_real32.h>
#include "ppolyval_spec.F90"

#include <numfort_real64.h>
#include "ppolyval_spec.F90"

    contains

#include <numfort_real32.h>
#include "ppolyval_impl.F90"

#include <numfort_real64.h>
#include "ppolyval_impl.F90"

end module