

#include <numfort.h>

module numfort_polynomial_ppoly_tensor

    use, intrinsic :: iso_fortran_env

    use numfort_core
    use numfort_common
    use numfort_common_status
    use numfort_common_input_checks
    use numfort_interpolate_common
    use numfort_interpolate_search

    implicit none
    private

    public :: ppoly2d
    public :: ppoly_get_ncoefs
    public :: ppoly_get_nknots
    public :: ppolyfit
    public :: ppolyval


    type :: ppoly_tensor_abc
        private
        integer :: degree = 0
    end type

    type, extends(ppoly_tensor_abc) :: ppoly2d
        private
        integer, dimension(2) :: nknots
    end type

    interface ppoly_get_ncoefs
        procedure ppoly2d_get_ncoefs
    end interface

    interface ppoly_get_nknots
        procedure ppoly2d_get_nknots
    end interface


#include <numfort_real32.h>
#include "ppoly_tensor_common_spec.F90"
#include "ppoly_tensor_fit_spec.F90"
#include "ppoly_tensor_val_spec.F90"

#include <numfort_real64.h>
#include "ppoly_tensor_common_spec.F90"
#include "ppoly_tensor_fit_spec.F90"
#include "ppoly_tensor_val_spec.F90"

    contains


pure function ppoly2d_get_ncoefs (self, n, k) result(res)
    type (ppoly2d), intent(in) :: self
    integer, intent(in), dimension(:), optional :: n
    integer, intent(in), optional :: k
    integer :: res

    integer :: lk, ln(2)

    lk = self%degree
    ln = self%nknots

    if (present(k)) lk = k
    if (present(n)) ln = n

    res = (ln(1)-1)*(ln(2)-1) * (lk+1) ** 2
end function


pure function ppoly2d_get_nknots (self, n, k) result(res)
    type (ppoly2d), intent(in) :: self
    integer, intent(in), dimension(:), optional :: n
    integer, intent(in), optional :: k
    integer :: res

    integer, dimension(2) :: ln

    ln = self%nknots
    if (present(n)) ln = n

    res = sum(ln)
end function

#include <numfort_real32.h>
#include "ppoly_tensor_common_impl.F90"
#include "ppoly_tensor_fit_impl.F90"
#include "ppoly_tensor_val_impl.F90"

#include <numfort_real64.h>
#include "ppoly_tensor_common_impl.F90"
#include "ppoly_tensor_fit_impl.F90"
#include "ppoly_tensor_val_impl.F90"

end module