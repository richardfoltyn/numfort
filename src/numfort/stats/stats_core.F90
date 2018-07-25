
#include <numfort.h>

module numfort_stats_core

    use, intrinsic :: iso_fortran_env
    use numfort_common

    implicit none
    private

    public :: mean
    public :: std
    public :: normalize
    public :: percentile


#include <numfort_real32.h>
#include "stats_core_spec.F90"

#include <numfort_real64.h>
#include "stats_core_spec.F90"


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

#include <numfort_real32.h>
#include "stats_core_impl.F90"

#include <numfort_real64.h>
#include "stats_core_impl.F90"

end module
