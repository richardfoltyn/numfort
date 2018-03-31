
#include <numfort.h>

module numfort_optimize_simplex
    !*  Wrapper module for CSIRO simplex implementation

    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_common_input_checks
    use numfort_common_workspace
    use numfort_optimize_result
    use numfort_optimize_interfaces
    use numfort_optimize_fwrapper

    use simplex_csiro, only: minim

    implicit none
    private

    public :: minimize_simplex

#include <numfort_real32.h>
#include "minimize_simplex_spec.F90"

#include <numfort_real64.h>
#include "minimize_simplex_spec.F90"


    contains

pure function map_iprint (i) result(k)
    !*  MAP_IPRINT converts NUMFORT print flags to values expected by
    !   wrapped routine.
    integer (NF_ENUM_KIND), intent(in) :: i
    integer :: k

    select case (i)
    case (NF_PRINT_MINIMAL)
        k = 0
    case (NF_PRINT_VERBOSE)
        k = 1
    case (NF_PRINT_ALL)
        k = 1
    case default
        ! do not print anything by default
        k = -1
    end select
end function


pure subroutine map_ifault (ifault, status, msg)
    !*  MAP_IFAULT maps the status code returned by the underlying implementation
    !   into corresponding NUMFORT status codes.
    integer, intent(in) :: ifault
        !!  Status code as returned by simplex routine
    type (status_t), intent(out) :: status
        !!  On exit, contains the corresponding NF status code
    character (len=*), intent(out), optional :: msg

    status = NF_STATUS_UNKNOWN
    if (ifault == 0) then
        status = NF_STATUS_OK
        if (present(msg)) msg = "Simplex: successful termination"
    else if (ifault == 1) then
        status = NF_STATUS_MAX_EVAL
        if (present(msg)) msg = "Simplex: max. number of function evaluations exceeded"
    else if (ifault == 3 .or. ifault == 4) then
        status = NF_STATUS_INVALID_ARG
        if (present(msg)) msg = "Simplex: invalid input argument(s)"
    end if
end subroutine


#include <numfort_real32.h>
#include "minimize_simplex_impl.F90"

#include <numfort_real64.h>
#include "minimize_simplex_impl.F90"


end module
