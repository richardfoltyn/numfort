
#include <numfort.h>


module numfort_polynomial_polyfit

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: vander
    use numfort_common
    use numfort_common_workspace
    use lapack_interfaces, only: GESV

    implicit none
    private

    public :: polyfit
    public :: polyfit_deriv

    interface polyfit
        procedure polyfit_real32, polyfit_1d_real32, &
            polyfit_real64, polyfit_1d_real64
    end interface
    
    interface polyfit_deriv
        procedure polyfit_deriv_1d_real32, polyfit_deriv_1d_real64
    end interface

    interface polyfit_exact
        procedure polyfit_exact_real32, polyfit_exact_real64
    end interface

    interface polyfit_check_input
        procedure polyfit_check_input_real32, polyfit_check_input_real64
    end interface

    contains


#include <numfort_real32.h>
#include "polyfit_impl.F90"

#include <numfort_real64.h>
#include "polyfit_impl.F90"

end module
