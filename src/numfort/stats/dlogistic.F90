#include "numfort.h"

module numfort_stats_dlogistic

    use, intrinsic :: iso_fortran_env

    implicit none
    private

    public :: pdf, cdf, ppf, rvs
    public :: logistic
   
#include "numfort_real32.h"
#include "cont/dlogistic_types.F90"

#include "numfort_real64.h"
#include "cont/dlogistic_types.F90"

    type (dlogistic_real64), protected :: logistic = dlogistic_real64(loc=0.0d0, scale=1.0d0)
        !*  Distribution object for the logistic distribution with default 
        !   location and scale parameters (double precision)

contains
   
#include "numfort_real32.h"
#include "cont/dlogistic_impl.F90"

#include "numfort_real64.h"
#include "cont/dlogistic_impl.F90"


end module
