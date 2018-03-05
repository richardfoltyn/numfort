
#include <numfort.h>

module numfort_common_testing
    !*  Module NUMFORT_COMMON_TESTING contains routines to perform various
    !   tests, eg. if two values are approximately close.

    use, intrinsic :: iso_fortran_env

    use numfort_common_shape

    implicit none

    private

    public :: all_close
    public :: is_close

#include <numfort_real32.h>
#include "include/testing_spec.F90"

#include <numfort_real64.h>
#include "include/testing_spec.F90"

    contains

#include <numfort_real32.h>
#include "include/testing_impl.F90"

#include <numfort_real64.h>
#include "include/testing_impl.F90"


end module
