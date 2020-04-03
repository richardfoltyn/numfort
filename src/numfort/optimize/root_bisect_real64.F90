

module numfort_optimize_bisect_real64

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_core, only: signum

    use numfort_optimize_result_real64
    use numfort_optimize_interfaces_real64
    use numfort_optimize_fwrapper_real64


    implicit none

    integer, parameter :: PREC = real64

    private

    public :: root_bisect

    integer, parameter :: MSG_LEN = 100

    interface root_bisect
        procedure root_bisect, root_bisect_args
    end interface

    contains

#include "root_bisect_impl.F90"

end module
