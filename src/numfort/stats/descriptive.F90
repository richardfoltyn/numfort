
#include "numfort.h"

module numfort_stats_desc

    use iso_fortran_env
    use numfort_common

    implicit none
    private

    interface mean
        module procedure mean_1d_real64, mean_2d_real64, &
            mean_1d_real32, mean_2d_real32
    end interface

    interface std
        module procedure std_1d_real64, std_2d_real64, &
            std_1d_real32, std_2d_real32
    end interface

    interface mean_std_check_input
        module procedure mean_std_check_input_real32, mean_std_check_input_real64
    end interface


    public :: mean, std
contains

! MEAN_STD_INIT sets the number of variables and number of observations
! of for a given 2-dimensional array.
pure subroutine mean_std_init (shp, dim, ldim, nvars, nobs, status)

    ! Input array shape; currently only 1d arrays with two elements
    ! are supported
    integer, intent(in), dimension(:) :: shp
    ! Dimension along with reduction should be performed (mean or std)
    integer, intent(in), optional :: dim
    ! On exit, contains dimension along which reduction should be applied
    ! (ldim), and number of variables (nvars), the number of observations (nobs)
    ! as well as an error code (status).
    integer, intent(out) :: ldim, status, nvars, nobs

    ldim = 1
    status = 0
    if (present(dim)) then
        if (dim /= 1 .and. dim /= 2) then
            status = 1
            return
        end if
        ldim = dim
    end if

    nvars = shp(3-ldim)
    nobs = shp(ldim)
end subroutine

#define __PREC real64
#include "descriptive_impl.F90"

#define __PREC real32
#include "descriptive_impl.F90"

end module
