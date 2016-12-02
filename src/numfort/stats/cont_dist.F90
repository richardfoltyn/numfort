
#include "dist_pdt.h"


module numfort_stats_cont_dist

    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    private


    type, abstract :: cont_dist __PDT_PARAM_DECL (PREC)
        __PDT_KIND_DECL(PREC, __DEFAULT_REAL_KIND)
    contains
        ! PDF method
        procedure (cont_dist_func_real64), deferred, private, pass :: pdf_real64
        generic, public :: pdf => pdf_real64

        ! CDF method
        procedure (cont_dist_func_real64), deferred, private, pass :: cdf_real64
        generic, public :: cdf => cdf_real64

        ! Random sample generation method
        procedure (cont_dist_real64), deferred, private, pass :: rvs_real64
        generic, public :: rvs => rvs_real64

    end type

    interface

#include "numfort_real64.h"
#include "cont/interfaces.F90"

#include "numfort_real32.h"
#include "cont/interfaces.F90"
    end interface

    public :: cont_dist

end module
