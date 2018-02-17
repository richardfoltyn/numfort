
#include <numfort.h>


module numfort_polynomial_polyroots

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_workspace

    implicit none
    private

    public :: polyroots

#include <numfort_real32.h>
#include "polyroots_spec.F90"

#include <numfort_real64.h>
#include "polyroots_spec.F90"


    contains

#include <numfort_real32.h>
#include "polyroots_impl.F90"

#include <numfort_real64.h>
#include "polyroots_impl.F90"

end module
