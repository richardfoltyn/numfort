
#include "dist_pdt.h"

module numfort_stats_dnorm

    use, intrinsic :: iso_fortran_env

    ! actual implementations from external libraries
    use cdf_normal_mod, only: cdf_normal
    use random, only: random_normal

    use numfort_core, only : PI
    use numfort_stats_cont_dist, only: cont_dist

    implicit none
    private

    real (real64), parameter :: NORM_CONST = 1/sqrt(2 * PI)

    !>  Implementation of the univariate normal distribution
    type, extends(cont_dist) :: dnorm
        real (__PDT_REAL_KIND) :: mean = 0.0
            !!  Mean of normal distribution
        real (__PDT_REAL_KIND) :: sd = 1.0
            !!  Standard deviation of normal distribution
    contains
        procedure, private, pass :: pdf_real64
        procedure, private, pass :: pdf_params_real64
        generic, public :: pdf => pdf_params_real64

        procedure, private, pass :: cdf_real64
        procedure, private, pass :: cdf_params_real64
        generic, public :: cdf => cdf_params_real64

        procedure, private, pass :: rvs_real64
        procedure, private, pass :: rvs_params_real64
        generic, public :: rvs => rvs_params_real64

        procedure, private, pass :: get_params_real64
        generic, private :: get_params => get_params_real64
    end type

    public :: norm, dnorm

    !>  Instance of the standard normal distribution
    type (dnorm), protected :: norm = dnorm(mean=0.0d0, sd=1.0d0)

contains

#include "numfort_real64.h"
#include "norm/routines.F90"


end module
