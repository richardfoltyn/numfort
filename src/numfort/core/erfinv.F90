

#include <numfort.h>

module numfort_core_erfinv
    !*  Module implements inverse of ERF() and some other related functions.
    !   The code here is a port of the C implementation in Scipy, which
    !   is an adapted version of the code in the Cephes math library.

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_core_constants


    implicit none
    private

    public :: erfinv
    public :: ndtri

#include <numfort_real32.h>
#include "erfinv_spec.F90"

#include <numfort_real64.h>
#include "erfinv_spec.F90"


    contains


#include <numfort_real32.h>
#include "erfinv_impl.F90"

#include <numfort_real64.h>
#include "erfinv_impl.F90"




end module
