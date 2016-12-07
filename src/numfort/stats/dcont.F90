
#include "numfort.h"
#include "dist_pdt.h"


module nf_stats_dcont

    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    private


    type, abstract :: dcont __PDT_PARAM_DECL (PREC)
        __PDT_KIND_DECL(PREC, __DEFAULT_REAL_KIND)
    contains
        ! PDF method
        procedure (x_fx_real64), deferred, private, pass :: pdf_real64
        generic, public :: pdf => pdf_real64

        ! CDF method
        procedure (x_fx_real64), deferred, private, pass :: cdf_real64
        generic, public :: cdf => cdf_real64

        ! Random sample generation method
        procedure (x_real64), deferred, private, pass :: rvs_real64
        generic, public :: rvs => rvs_real64

    end type

    interface

#define __PREC real64
#include "cont/interfaces.F90"

#define __PREC real32
#include "cont/interfaces.F90"
    end interface

    public :: dcont

end module
