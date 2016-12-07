
#include <numfort.h>
#include "dist_pdt.h"

module nf_stats_ddisc

    use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
    implicit none
    private


    type, abstract :: ddisc __PDT_PARAM_DECL_BOTH (PREC, INTSIZE)
        __PDT_KIND_DECL(PREC, __DEFAULT_REAL_KIND)
        __PDT_KIND_DECL(INTSIZE, __DEFAULT_INT_KIND)
    contains
        ! pmf method
        procedure (__APPEND2(x_fx,__DEFAULT_INT_KIND,real64)), deferred, private, pass :: &
            __APPEND2(pmf,__DEFAULT_INT_KIND,real64)
        generic, public :: pmf => __APPEND2(pmf,__DEFAULT_INT_KIND,real64)

        ! CDF method
        ! procedure, private, pass :: cdf_scalar
        ! procedure (disc_dist_func_array), deferred, private, pass :: cdf_array
        ! generic, public :: cdf => cdf_scalar, cdf_array

        ! Random sample generation method
        procedure (__APPEND2(x,__DEFAULT_INT_KIND,real64)), deferred, private, pass :: &
            __APPEND2(rvs,__DEFAULT_INT_KIND,real64)
        generic, public :: rvs => __APPEND2(rvs,__DEFAULT_INT_KIND,real64)

    end type

    interface
#define __PREC real64
#define __INTSIZE __DEFAULT_INT_KIND
#include "disc/interfaces.F90"
    end interface

    public :: ddisc

end module
