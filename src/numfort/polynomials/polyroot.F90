
#include <numfort.h>


module numfort_polynomial_polyroot

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_workspace

    implicit none
    private

    public :: polyroot

    interface polyroot
        procedure polyroot_real32, polyroot_real64
    end interface

    interface polyroot_quad
        procedure polyroot_quad_real32, polyroot_quad_real64
    end interface

    interface polyroot_check_input
        procedure polyroot_check_input_real32, polyroot_check_input_real64
    end interface

    contains

#include <numfort_real32.h>
#include "polyroot_impl.F90"

#include <numfort_real64.h>
#include "polyroot_impl.F90"

end module
