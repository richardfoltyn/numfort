
#include "numfort.h"

module numfort_optimize_brent

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_optimize_result

    implicit none
    private

    public :: root_brentq

    abstract interface
        subroutine fcn_real32 (x, fx)
            import real32
            real (real32), intent(in) :: x
            real (real32), intent(out) :: fx
        end subroutine

        subroutine fcn_real64 (x, fx)
            import real64
            real (real64), intent(in) :: x
            real (real64), intent(out) :: fx
        end subroutine
    end interface

    interface root_brentq
        module procedure root_brentq_real32, root_brentq_real64
    end interface

    interface root_brentq_impl
        module procedure root_brentq_impl_real32, root_brentq_impl_real64
    end interface

contains

#define __PREC real32
#include "root_brent_impl.F90"
#undef __PREC

#define __PREC real64
#include "root_brent_impl.F90"
#undef __PREC

end module
