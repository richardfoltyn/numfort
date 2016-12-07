
#include <numfort.h>
#include "dist_pdt.h"

module nf_stats_drandint

    use iso_fortran_env
    use nf_stats_ddisc

    implicit none
    private


    type, extends(ddisc) :: drandint
        integer (__PDT_INT_KIND) :: low = 0
        integer (__PDT_INT_KIND) :: high = 1
    contains
        procedure, private, pass :: __APPEND2(pmf,__DEFAULT_INT_KIND,real64)
        procedure, private, pass :: __APPEND2(pmf_params,__DEFAULT_INT_KIND,real64)
        generic, public :: pmf => __APPEND2(pmf_params,__DEFAULT_INT_KIND,real64)

        procedure, private, pass :: __APPEND2(rvs,__DEFAULT_INT_KIND,real64)
        procedure, private, pass :: __APPEND2(rvs_params,__DEFAULT_INT_KIND,real64)
        generic, public :: rvs => __APPEND2(rvs_params,__DEFAULT_INT_KIND,real64)
    end type

    type (drandint), parameter :: randint = drandint(low=0, high=1)

    public :: randint, drandint

contains

#define __PREC real64
#define __INTSIZE __DEFAULT_INT_KIND
#include "randint/routines.F90"

end module
