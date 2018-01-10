#include "numfort.h"

module numfort_stats_dnorm

    use, intrinsic :: iso_fortran_env

    ! actual implementations from external libraries
    use cdf_normal_mod, only: cdf_normal
    use random, only: random_normal

    use numfort_core, only : PI

    implicit none
    private

    public :: pdf, cdf, rvs
    public :: norm, dnorm_real64

    real (real64), parameter :: NORM_CONST = 1/sqrt(2 * PI)

#include "numfort_real64.h"
#include "cont/dnorm_types.F90"

    type (dnorm_real64), parameter :: norm = dnorm_real64(loc=0.0d0, scale=1.0d0)
        !*  Instance of the standard normal distribution (double precision)

contains

#include "numfort_real64.h"
#include "cont/dnorm_impl.F90"

end module
