
#include "numfort.h"

module numfort_stats_dgenpareto
    !*  Implementation of Generalized Pareto Distribution
    !   with CDF
    !       F(x) = 1 - (1 + xi(x-mu)/sigma) ** (-1/xi)
    !   where mu, sigma and xi are the location, scale
    !   and shape parameters (as listed in Wikipedia
    !   entry for GPD).

    use iso_fortran_env
    implicit none

    private

    public :: cdf, pdf


#include "numfort_real32.h"
#include "cont/dgenpareto_types.F90"

#include "numfort_real64.h"
#include "cont/dgenpareto_types.F90"


contains

#include "numfort_real32.h"
#include "cont/dgenpareto_impl.F90"

#include "numfort_real64.h"
#include "cont/dgenpareto_impl.F90"

end module
