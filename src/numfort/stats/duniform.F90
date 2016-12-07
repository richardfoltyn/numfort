#include "numfort.h"
#include "dist_pdt.h"

module nf_stats_duniform

    use, intrinsic :: iso_fortran_env
    use nf_stats_dcont, only: dcont

    implicit none
    private

    !>  Implementation of the continuous univariate uniform distribution
    !   on an interval [low,high]
    type, extends(dcont) :: duniform
        real (__PDT_REAL_KIND) :: low = 0.0
            !!  Lower bound of support
        real (__PDT_REAL_KIND) :: high = 1.0
            !!  Upper bound of support
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
        generic :: get_params => get_params_real64
    end type

    interface pdf_impl
        module procedure pdf_impl_real64
    end interface

    interface cdf_impl
        module procedure cdf_impl_real64
    end interface

    public :: uniform, duniform

    !>  Instance of the standard uniform distribution on [0,1].
    type (duniform), protected :: uniform = duniform(low=0.0d0, high=1.0d0)

contains

#define __PREC real64
#include "uniform/routines.F90"

end module
