


#include <numfort.h>

module numfort_polynomial_ppoly
    !*  Module implements fitting, evaluation and interpolation of
    !   piecewise polynomials.

    use, intrinsic :: iso_fortran_env

    use numfort_core
    use numfort_common
    use numfort_common_status
    use numfort_common_input_checks
    use numfort_interpolate

    implicit none
    private

    public :: ppoly
    public :: ppoly_bernstein
    public :: ppoly_get_nknots
    public :: ppoly_get_ncoefs
    public :: ppolyfit
    public :: ppolyval
    public :: ppolyder

    ! The following routines contain the implementations of the more
    ! userfriendly routines exported above, but do not perform any input checks.
    ! For NUMFORT-internal use only!
    public :: ppolyval_impl
    public :: bernstein_fit_deriv_impl
    public :: ppolyder_impl

    type :: ppoly_abc
        private
        integer :: degree = 0
            !*  Polynomial degree of each piecewise polynomial
        integer :: nknots = 0
            !*  Number of knots (breakpoints)
    end type

    type, extends(ppoly_abc) :: ppoly
        !*  Piecewise polynomial, default basis (ie. power basis)
    end type

    type, extends(ppoly_abc) :: ppoly_bernstein
        !*  Bernstein polynomial basis
    end type

    interface ppoly_get_nknots
        procedure power_get_nknots, bernstein_get_nknots
    end interface

    interface ppoly_get_ncoefs
        procedure power_get_ncoefs, bernstein_get_ncoefs
    end interface

#include <numfort_real32.h>
#include "ppoly_common_spec.F90"
#include "ppolyfit_spec.F90"
#include "ppolyval_spec.F90"
#include "ppolyder_spec.F90"

#include <numfort_real64.h>
#include "ppoly_common_spec.F90"
#include "ppolyfit_spec.F90"
#include "ppolyval_spec.F90"
#include "ppolyder_spec.F90"

    contains

pure function power_get_ncoefs (self, n, k) result(res)
    !*  POWER_GET_NCOEFS returns the required array size to
    !   store the data of a piecewise polynomial of given degree K with N break
    !   points using the power basis.
    type (ppoly), intent(in) :: self
    integer, intent(in), optional :: n
    integer, intent(in), optional :: k
    integer :: res

    integer :: lk, ln

    lk = self%degree
    ln = self%nknots
    if (present(k)) lk = k
    if (present(n)) ln = n

    res = (ln-1) * (lk + 1)
end function


pure function power_get_nknots (self, n, k) result(res)
    type (ppoly), intent(in) :: self
    integer, intent(in), optional :: n
    integer, intent(in), optional :: k
    integer :: res

    integer :: ln

    ln = self%nknots
    if (present(n)) ln = n
    res = ln
end function


pure function bernstein_get_ncoefs (self, n, k) result(res)
    type (ppoly_bernstein), intent(in) :: self
    integer, intent(in), optional :: n
    integer, intent(in), optional :: k
    integer :: res

    integer :: lk, ln

    lk = self%degree
    ln = self%nknots
    if (present(k)) lk = k
    if (present(n)) ln = n
    res = (ln-1) * (lk + 1)
end function


pure function bernstein_get_nknots (self, n, k) result(res)
    type (ppoly_bernstein), intent(in) :: self
    integer, intent(in), optional :: n
    integer, intent(in), optional :: k
    integer :: res

    integer :: ln

    ln = self%nknots
    if (present(n)) ln = n
    res = ln
end function


#include <numfort_real32.h>
#include "ppoly_common_impl.F90"
#include "ppolyfit_impl.F90"
#include "ppolyval_impl.F90"
#include "ppolyder_impl.F90"

#include <numfort_real64.h>
#include "ppoly_common_impl.F90"
#include "ppolyfit_impl.F90"
#include "ppolyval_impl.F90"
#include "ppolyder_impl.F90"

end module