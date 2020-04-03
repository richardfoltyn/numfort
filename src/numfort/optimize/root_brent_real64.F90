

module numfort_optimize_brent_real64

    use, intrinsic :: iso_fortran_env

    use numfort_common_status
    use numfort_optimize_result_real64
    use numfort_optimize_interfaces_real64

    implicit none

    integer, parameter :: PREC = real64

    private

    public :: root_brentq

    interface root_brentq
        procedure root_brentq
    end interface


    contains

#include "root_brent_impl.F90"

end module
