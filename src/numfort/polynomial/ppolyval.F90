
#include <numfort.h>

module numfort_polynomial_ppolyval

    use, intrinsic :: iso_fortran_env

    use numfort_core
    use numfort_common
    use numfort_common_status
    use numfort_common_workspace
    use numfort_interpolate_common

    use numfort_common
    use numfort_core

    implicit none
    private

    public :: ppolyval
    public :: ppoly_get_ncoef

#include <numfort_real32.h>
#include "ppolyval_spec.F90"

#include <numfort_real64.h>
#include "ppolyval_spec.F90"

    contains

pure function ppoly_get_ncoef (n, k) result(res)
    !*  PPOLY_GET_NCOEF returns the required minimum array size to
    !   store the data of a piecewise polynomial of given degree K with N break
    !   points.
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer :: res

    res = (n-1) * (k + 1)
end function


#include <numfort_real32.h>
#include "ppolyval_impl.F90"

#include <numfort_real64.h>
#include "ppolyval_impl.F90"

end module