
#include "numfort.h"

module numfort_integrate_quadpack_real64
    
    use, intrinsic :: iso_fortran_env
    use quadpack
    use numfort_common_status
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_integrate_quadpack_common

    implicit none

    private

    public :: quad

    abstract interface
        function f_integrand (x) result(fx)
            import real64
            real (real64), intent(in) :: x
            real (real64) :: fx
        end function
    end interface

    interface quad
        procedure quad
    end interface

    integer, parameter :: PREC = real64

contains

#include "quad_impl.F90"

end module
