
#include <numfort.h>

module numfort_polynomial_polyval

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: vander
    use numfort_common
    use numfort_common_workspace

    implicit none
    private


    contains

#include <numfort_real32.h>
#include "polyval_impl.F90"

end module
