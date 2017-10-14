
#include <numfort.h>

module numfort_optimize_interfaces

    use, intrinsic :: iso_fortran_env
    implicit none

    abstract interface
#include <numfort_real32.h>
#include "interfaces_spec.F90"

#include <numfort_real64.h>
#include "interfaces_spec.F90"

    end interface

end module
