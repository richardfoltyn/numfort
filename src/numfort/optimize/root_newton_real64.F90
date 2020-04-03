

module numfort_optimize_newton_real64

    use, intrinsic :: iso_fortran_env

    use numfort_common_status
    use numfort_common_swap
    use numfort_common_input_checks
    use numfort_core, only: signum

    use numfort_optimize_result_real64
    use numfort_optimize_interfaces_real64
    use numfort_optimize_fwrapper_real64

    implicit none

    integer, parameter :: PREC = real64

    private

    public :: root_newton
    public :: root_halley
    public :: root_newton_bisect

    abstract interface
        subroutine fcn_der2_args (x, args, fx, fpx, fppx)
            import PREC
            real (PREC), intent(in) :: x
            real (PREC), intent(inout), dimension(:) :: args
            real (PREC), intent(out) :: fx, fpx, fppx
        end subroutine

        subroutine fcn_der2 (x, fx, fpx, fppx)
            import PREC
            real (PREC), intent(in) :: x
            real (PREC), intent(out) :: fx, fpx, fppx
        end subroutine
    end interface

    interface root_newton
        procedure root_newton, root_newton_args, &
            root_newton_jac, root_newton_jac_args, &
            root_newton_fcn_jac, root_newton_fcn_jac_args
    end interface

    interface root_halley
        procedure root_halley, root_halley_args
    end interface

    interface root_newton_bisect
        procedure newton_bisect, newton_bisect_args, &
            newton_bisect_jac, newton_bisect_jac_args, &
            newton_bisect_fcn_jac, newton_bisect_fcn_jac_args
    end interface


    integer, parameter :: MSG_LEN = 100

    contains

#include "root_newton_impl.F90"

end module
