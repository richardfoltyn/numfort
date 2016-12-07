
! ------------------------------------------------------------------------------
! GET_PARAMS method
pure subroutine __APPEND(get_params,__PREC) (self, mean, sd, mean_out, sd_out)
    integer, parameter :: PREC = __PREC
    class (dnorm __PDT_PARAM_DECL(PREC)), intent(in) :: self
    real (PREC), intent(in), optional :: mean, sd
    real (PREC), intent(out) :: mean_out, sd_out

    if (present(mean)) then
        mean_out = mean
    else
        mean_out = real(self%mean, PREC)
    end if

    if (present(sd)) then
        sd_out = sd
    else
        sd_out = real(self%sd, PREC)
    end if
end subroutine

! ------------------------------------------------------------------------------
! PDF method
impure elemental subroutine __APPEND(pdf_impl,__PREC) (x, fx, mean, sd)
    integer, parameter :: PREC = __PREC
    ! class (dnorm __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"
#include "norm/params.F90"
    intent (out) :: fx

    real (PREC) :: b

    b = 2 * (sd ** 2)
    fx = NORM_CONST/sd * exp(-((x-mean) ** 2) / b)
end subroutine

impure elemental function __APPEND(pdf_params,__PREC) (self, x, mean, sd) result(fx)
    integer, parameter :: PREC = __PREC
    class (dnorm __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"
#include "norm/params.F90"

    call pdf_impl (x, fx, mean, sd)
end function

impure elemental function __APPEND(pdf,__PREC) (self, x) result(fx)
    integer, parameter :: PREC = __PREC
    class (dnorm __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"

    call pdf_impl (x, fx, mean=self%mean, sd=self%sd)
end function

! ------------------------------------------------------------------------------
! CDF method

impure elemental subroutine __APPEND(cdf_impl,__PREC) (x, fx, mean, sd)
    integer, parameter :: PREC = __PREC
#include "cont/spec_func.F90"
#include "norm/params.F90"
    intent (out) :: fx

    ! CDFLIB90 argument
    integer, parameter :: which = 1
    ! CDFLIB90 refuses to operate on domain outside of [-immense, immense]
    real (PREC), parameter :: immense = 1.0e100_PREC

    ! Check whether we are outside of supported domain
    if ((x-mean) < -immense) then
        fx = 0.0_PREC
        return
    else if ((x-mean) > immense) then
        fx = 1.0_PREC
        return
    end if

    ! call CDFLIB90 routine to actually do the work
    call cdf_normal (which, cum=fx, x=x, mean=mean, sd=sd)

end subroutine

impure elemental function __APPEND(cdf_params,__PREC) (self, x, mean, sd) result(fx)
    integer, parameter :: PREC = __PREC
    class (dnorm __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"
#include "norm/params.F90"

    call cdf_impl (x, fx, mean, sd)
end function

impure elemental function __APPEND(cdf,__PREC) (self, x) result(fx)
    integer, parameter :: PREC = __PREC
    class (dnorm __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec_func.F90"

    call cdf_impl (x, fx, mean=self%mean, sd=self%sd)
end function

! ------------------------------------------------------------------------------
! RVS method

impure elemental subroutine __APPEND(rvs_params,__PREC) (self, x, mean, sd)
    integer, parameter :: PREC = __PREC
    class (dnorm __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec.F90"
#include "norm/params.F90"

    real (PREC) :: z

    ! random draw from std. normal distribution
    z = random_normal ()
    x = mean + z*sd
end subroutine

impure elemental subroutine __APPEND(rvs,__PREC) (self, x)
    integer, parameter :: PREC = __PREC
    class (dnorm __PDT_PARAM_DECL(PREC)), intent(in) :: self
#include "cont/spec.F90"

    call self%rvs (x, mean=self%mean, sd=self%sd)
end subroutine
