
#include <numfort.h>

module numfort_interpolate_pchip
    !*  Module implements fitting and evaluation using the PCHIP algorithm
    !   (Piecewise Cubic Hermite Intepolating polynomials) that preserves
    !   the monotonicity of monotonic input data.
    !
    !   References:
    !       Fritsch, F. N. & Carlson, R. E. (1980):  Monotone Piecewise
    !           Cubic Interpolation. SIAM Journal on Numerical Analysis (17) 238-246
    use, intrinsic :: iso_fortran_env

    use numfort_core
    use numfort_common
    use numfort_common_status
    use numfort_common_workspace
    use numfort_interpolate_common
    use numfort_interpolate_search

    implicit none

    private

    public :: interp_pchip_get_ncoef
    public :: interp_pchip_fit
    public :: interp_pchip_eval

    integer, parameter :: POLY_DEGREE = 3
        !*  Degree of interpolating polynomial (fixed, cannot be changed
        !   by user).
    integer, parameter :: MAX_ORDER = 2
        !*  Max. order of derivative that can be evaluated

    integer, parameter :: IDX_Y = 1, IDX_D = 2, IDX_C = 3, IDX_B = 4
        !*   Positions of individual coefficients within each interval block

#include "numfort_real64.h"
#include "interp_pchip_spec.F90"



    contains


pure function interp_pchip_get_ncoef (n) result(res)
    !*  INTERP_PCHIP_GET_SIZE returns the required minimum array size to
    !   store PCHIP interpolation data for N break points.
    integer, intent(in) :: n
    integer :: res

    ! Need to store 4 data points per interval
    res = (n-1) * 4
end function

#include "numfort_real64.h"
#include "interp_pchip_impl.F90"


end module