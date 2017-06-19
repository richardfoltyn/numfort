
#include "numfort.h"

module numfort_integrate_quadpack
    
    use, intrinsic :: iso_fortran_env
    use quadpack
    use numfort_common_enums
    use numfort_common_status
    use numfort_common_workspace

    implicit none

    private

    public :: quad

    ! Custom error codes for QUADPACK wrappers
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_MAX_INTERVALS = ishft(1, 1)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_INTEGRAND_ERROR = ishft(1, 2)

    abstract interface
        function f_integrand_real32 (x) result(fx)
            import real32
            integer, parameter :: PREC = real32
            real (PREC), intent(in) :: x
            real (PREC) :: fx
        end function
        function f_integrand_real64 (x) result(fx)
            import real64
            integer, parameter :: PREC = real64
            real (PREC), intent(in) :: x
            real (PREC) :: fx
        end function
    end interface


    interface quad
        module procedure quad_real64
    end interface

    interface quad_check_input
        module procedure quad_check_input_real64
    end interface

contains

subroutine map_qagse_status (ier, status)
    !*  MAP_QAGSE_STATUS maps the integer value returned by 
    !   QUADPACK's QAGSE routine to NUMFORT status codes.

    integer, intent(in) :: ier
    type (status_t), intent(in out) :: status

    select case (ier)
    case (0)
        status = NF_STATUS_OK
    case (1)
        status = NF_STATUS_MAX_INTERVALS
    case (2)
        status = NF_STATUS_NOT_CONVERGED
    case (3)
        status = NF_STATUS_INTEGRAND_ERROR
    case (4)
        status = NF_STATUS_NOT_CONVERGED
    case (5)
        status = NF_STATUS_INTEGRAND_ERROR
    case (6)
        status = NF_STATUS_INVALID_ARG
    end select

    status%code_orig = ier
end subroutine

#include "numfort_real64.h"
#include "quad_impl.F90"

end module
