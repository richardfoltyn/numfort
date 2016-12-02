impure elemental subroutine __APPEND_PREC(cont_dist_func) (self, x, fx)
    import __PREC
    import cont_dist
    integer, parameter :: PREC = __PREC
    class (cont_dist __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"
end subroutine

impure elemental subroutine __APPEND_PREC(cont_dist) (self, x)
    import __PREC
    import cont_dist
    integer, parameter :: PREC = __PREC
    class (cont_dist __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec.F90"
end subroutine
