
! ------------------------------------------------------------------------------
! GET_PARAMS method
pure subroutine __APPEND_PREC(get_params) (self, low, high, low_out, high_out)
    integer, parameter :: PREC = __PREC
    class (duniform __PDT_PARAM_DECL(PREC)), intent(in) :: self
    real (PREC), intent(in), optional :: low, high
    real (PREC), intent(out) :: low_out, high_out

    if (present(low)) then
        low_out = low
    else
        low_out = real(self%low, PREC)
    end if

    if (present(high)) then
        high_out = high
    else
        high_out = real(self%high, PREC)
    end if
end subroutine

! ------------------------------------------------------------------------------
! PDF method

impure elemental subroutine __APPEND_PREC(pdf_impl) (x, fx, low, high)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
#include "uniform/params.F90"

    fx = 1 / (high-low)
    ! clip PDF value outside of valid range
    if (x < low) then
        fx = 0.0_PREC
    else if (x > high) then
        fx = 0.0_PREC
    end if
end subroutine

impure elemental function __APPEND_PREC(pdf_params) (self, x, low, high) result(fx)
    integer, parameter :: PREC = __PREC
    class (duniform __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"
#include "uniform/params.F90"

    call pdf_impl (x, fx, low, high)
end function

impure elemental function __APPEND_PREC(pdf) (self, x) result(fx)
    integer, parameter :: PREC = __PREC
    class (duniform __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"

    call pdf_impl (x, fx, low=self%low, high=self%high)
end function

! ------------------------------------------------------------------------------
! CDF method

impure elemental subroutine __APPEND_PREC(cdf_impl) (x, fx, low, high)
    integer, parameter :: PREC = __PREC
    ! class (duniform __PDT_PARAM_DECL(PREC)), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
#include "uniform/params.F90"

    fx = (x - low) / (high-low)
    ! adjust CDF values outside of valid range
    if (x < low) then
        fx = 0.0_PREC
    else if (x > high) then
        fx = 1.0_PREC
    end if
end subroutine

impure elemental function __APPEND_PREC(cdf_params) (self, x, low, high) result(fx)
    integer, parameter :: PREC = __PREC
    class (duniform __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"
#include "uniform/params.F90"

    call cdf_impl (x, fx, low, high)
end function

impure elemental function __APPEND_PREC(cdf) (self, x) result(fx)
    integer, parameter :: PREC = __PREC
    class (duniform __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"

    call cdf_impl (x, fx, low=self%low, high=self%high)
end function

! ------------------------------------------------------------------------------
! RVS method

impure elemental subroutine __APPEND_PREC(rvs_impl) (x, low, high)
    integer, parameter :: PREC = __PREC
    ! class (duniform __PDT_PARAM_DECL(PREC)), intent(in) :: self
    real (PREC), intent(out) :: x
#include "uniform/params.F90"

    call random_number (x)
    x = x * (high-low) + low
end subroutine

impure elemental function __APPEND_PREC(rvs_params) (self, low, high) result(x)
    integer, parameter :: PREC = __PREC
    class (duniform __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec.F90"
#include "uniform/params.F90"

    call rvs_impl (x, low, high)
end function

impure elemental function __APPEND_PREC(rvs) (self) result(x)
    integer, parameter :: PREC = __PREC
    class (duniform __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec.F90"

    call rvs_impl (x, low=self%low, high=self%high)
end function
