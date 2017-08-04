
#include "numfort.h"

module numfort_integrate_basic

    use, intrinsic :: iso_fortran_env
    use numfort_common_status

    implicit none

    private

    public :: trapezoid 

    interface trapezoid
        module procedure trapezoid_real32, trapezoid_real64, &
            trapezoid_array_real32, trapezoid_array_real64
    end interface

    interface check_inputs
        module procedure check_inputs_real32, check_inputs_real64, &
            check_inputs_array_real32, check_inputs_array_real64
    end interface

contains

#include "numfort_real32.h"
#include "trapezoid_impl.F90"

#include "numfort_real64.h"
#include "trapezoid_impl.F90"

end module

