
#include "numfort.h"

module numfort_stats_dfrechet
    !*  Implementation of the Frechet distribution with CDF
    !       F(x) = exp(-(x-m)/s)**(-alpha))
    !   where m, s and alpha are the location, scale
    !   and shape parameters (using the notation on Wikipedia's Frechet
    !   distribution page).

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env
    implicit none

    private

    public :: cdf, pdf, ppf, mean


#include "numfort_real32.h"
#include "cont/dfrechet_types.F90"

#include "numfort_real64.h"
#include "cont/dfrechet_types.F90"


contains

#include "numfort_real32.h"
#include "cont/dfrechet_impl.F90"

#include "numfort_real64.h"
#include "cont/dfrechet_impl.F90"

end module
