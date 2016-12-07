

! ------------------------------------------------------------------------------
! PMF method

impure elemental function __APPEND2(pmf,__INTSIZE,__PREC) (self, x) result(fx)
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = __INTSIZE
    class (drandint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
#include "disc/spec_func.F90"

end function

impure elemental function __APPEND2(pmf_params,__INTSIZE,__PREC) (self, x, low, high) result(fx)
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = __INTSIZE
    class (drandint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
#include "disc/spec_func.F90"
#include "randint/params.F90"

end function

! ------------------------------------------------------------------------------
! RVS method

impure elemental subroutine __APPEND2(rvs_params,__INTSIZE,__PREC) (self, x, low, high)
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = __INTSIZE
    class (drandint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
#include "disc/spec.F90"
#include "randint/params.F90"

end subroutine

impure elemental subroutine __APPEND2(rvs,__INTSIZE,__PREC) (self, x)
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = __INTSIZE
    class (drandint __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
#include "disc/spec.F90"

end subroutine
