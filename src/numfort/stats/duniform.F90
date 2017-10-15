#include "numfort.h"

module numfort_stats_duniform

    use, intrinsic :: iso_fortran_env

    implicit none
    private

    public :: pdf, cdf, rvs
   
#include "numfort_real32.h"
#include "cont/duniform_types.F90"

#include "numfort_real64.h"
#include "cont/duniform_types.F90"

    type (duniform_real64), protected, public :: &
        uniform = duniform_real64(loc=0.0d0, scale=1.0d0)
        !*  Distribution object for the continuous standard uniform
        !   distribution (double precision)

contains
   
#include "numfort_real32.h"
#include "cont/duniform_impl.F90"

#include "numfort_real64.h"
#include "cont/duniform_impl.F90"


end module
