impure elemental function __APPEND2(x_fx,__INTSIZE,__PREC) (self, x) result(fx)
    import __PREC, __INTSIZE
    import ddisc
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = __INTSIZE
    class (ddisc __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
#include "disc/spec_func.F90"
end function

impure elemental subroutine __APPEND2(x,__INTSIZE,__PREC) (self, x)
    import __PREC, __INTSIZE
    import ddisc
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = __INTSIZE
    class (ddisc __PDT_PARAM_DECL_BOTH(PREC, INTSIZE)), intent(in) :: self
#include "disc/spec.F90"
end subroutine
