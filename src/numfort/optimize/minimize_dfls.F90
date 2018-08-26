
#include <numfort.h>

module numfort_optimize_dfls

    use, intrinsic :: iso_fortran_env

    use numfort_common_enums
    use numfort_common_status
    use numfort_common_input_checks
    use numfort_common_workspace

    use numfort_optimize_interfaces
    use numfort_optimize_fwrapper
    use numfort_optimize_result

    use newuoa2_real32
    use newuoa2_real64

    implicit none

    private

    public :: minimize_dfls

#include <numfort_real32.h>
#include "minimize_dfls_spec.F90"

#include <numfort_real64.h>
#include "minimize_dfls_spec.F90"

    contains


pure function get_print_level (val) result(res)
    integer (NF_ENUM_KIND), intent(in) :: val
    integer :: res

    select case (val)
    case (NF_PRINT_ALL)
        res = 3
    case (NF_PRINT_VERBOSE)
        res = 2
    case (NF_PRINT_MINIMAL)
        res = 1
    case default
        res = 0
    end select
end function


#include <numfort_real32.h>
#include "minimize_dfls_impl.F90"

#include <numfort_real64.h>
#include "minimize_dfls_impl.F90"

end module