
#include <numfort.h>

module numfort_stats_core

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: argsort
    use numfort_common
    use numfort_interpolate

    implicit none
    private

    public :: mean
    public :: std
    public :: normalize
    public :: quantile


#include <numfort_real32.h>
#include "stats_core_spec.F90"

#include <numfort_real64.h>
#include "stats_core_spec.F90"

    integer (NF_ENUM_KIND), parameter :: NF_STATS_QUANTILE_LINEAR = 1
    integer (NF_ENUM_KIND), parameter :: NF_STATS_QUANTILE_LOWER = 2
    integer (NF_ENUM_KIND), parameter :: NF_STATS_QUANTILE_HIGHER = 4
    integer (NF_ENUM_KIND), parameter :: NF_STATS_QUANTILE_NEAREST = 8
    integer (NF_ENUM_KIND), parameter :: NF_STATS_QUANTILE_MIDPOINT = 16

    contains

pure subroutine mean_std_init (shp, dim, ldim, nvars, nobs, status)
    !*  MEAN_STD_INIT sets the number of variables and number of observations
    !   of for a given 2-dimensional array.
    integer, intent(in), dimension(:) :: shp
        !*  Input array shape; currently only 1d arrays with two elements
        !   are supported
    integer, intent(in), optional :: dim
        !*  Optional user-provided dimension along with reduction should be
        !   performed (mean or std)
    integer, intent(out) :: ldim
        !*  Dimension along with reduction should be applied, based on
        !   user-provided arguments or default rules.
    integer, intent(out) :: nvars
        !*  Number of variables
    integer, intent(out) :: nobs
        !*  Number of observations
    type (status_t), intent(out) :: status
        !*  Status flag (either NF_STATUS_OK or NF_STATUS_INVALID_ARG)

    ldim = 1
    status = NF_STATUS_OK
    if (present(dim)) then
        if (dim /= 1 .and. dim /= 2) then
            status = NF_STATUS_INVALID_ARG
            return
        end if
        ldim = dim
    end if

    nvars = shp(3-ldim)
    nobs = shp(ldim)
end subroutine


pure function quantile_interp_to_enum (interp) result(res)
    !*  PERCENTILE_INTEPR_TO_ENUM converts the character-type interpolation
    !   method to its integer representation.
    character (*), intent(in) :: interp
    integer (NF_ENUM_KIND) :: res

    character (10) :: linterp

    linterp = interp
    call lower (linterp)

    select case (linterp)
    case ("linear")
        res = NF_STATS_QUANTILE_LINEAR
    case ("lower")
        res = NF_STATS_QUANTILE_LOWER
    case ("higher")
        res = NF_STATS_QUANTILE_HIGHER
    case ("nearest")
        res = NF_STATS_QUANTILE_NEAREST
    case ("midpoint")
        res = NF_STATS_QUANTILE_MIDPOINT
    case default
        res = 0
    end select
end function



#include <numfort_real32.h>
#include "stats_core_impl.F90"

#include <numfort_real64.h>
#include "stats_core_impl.F90"

end module
