impure elemental function __APPEND(x_fx,__PREC) (self, x) result(fx)
    import __PREC
    import dcont
    integer, parameter :: PREC = __PREC
    class (dcont __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"
end function

impure elemental subroutine __APPEND(x,__PREC) (self, x)
    import __PREC
    import dcont
    integer, parameter :: PREC = __PREC
    class (dcont __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec.F90"
end subroutine
