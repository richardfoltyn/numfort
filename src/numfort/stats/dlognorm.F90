

#include <numfort.h>

module numfort_stats_dlognorm
    !*  Module implements various routines for the log-normal distribution.
    !
    !   Note that the distribution is parametrized in terms of the location
    !   and scale parameters of the underlying normal distribution, even
    !   though these are not location or scale parameters of the log-normal
    !   distribution per se (thus the implementation deviates from Scipy's
    !   lognorm class).

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    ! actual implementations from external libraries
    use random, only: random_normal, random_normal_real32, random_normal_real64

    use numfort_core

    implicit none
    private

    public :: pdf
    public :: cdf
    public :: ppf

#include <numfort_real64.h>
#include "cont/dlognorm_spec.F90"

    contains

#include <numfort_real64.h>
#include "cont/dlognorm_impl.F90"

end module
