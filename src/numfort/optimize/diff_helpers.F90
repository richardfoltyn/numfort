
#include <numfort.h>

module numfort_optimize_diff
    !*  Module contains helper routines to numerically differentiate
    !   (multivariate) functions.

    use, intrinsic :: iso_fortran_env

    use numfort_optimize_interfaces

    implicit none

    private
    public :: num_diff

#include <numfort_real32.h>
#include "diff_helpers_spec.F90"

#include <numfort_real64.h>
#include "diff_helpers_spec.F90"


    contains

#include <numfort_real32.h>
#include "diff_helpers_impl.F90"

#include <numfort_real64.h>
#include "diff_helpers_impl.F90"

end module
