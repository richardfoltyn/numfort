

#include <numfort.h>

module numfort_core_cephes_stats
    !*  Module implements various special functions useful for statistics,
    !   such as the normal CDF, inverse CDF and the inverse ERF function.
    !   The code here is a port of the C implementation in Scipy, which
    !   is an adapted version of the code in the Cephes math library.

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_core_constants


    implicit none
    private

    public :: erfinv
    public :: ndtr
    public :: ndtri

#include <numfort_real32.h>
#include "cephes_stats_spec.F90"

#include <numfort_real64.h>
#include "cephes_stats_spec.F90"


    contains


#include <numfort_real32.h>
#include "cephes_stats_impl.F90"

#include <numfort_real64.h>
#include "cephes_stats_impl.F90"




end module
