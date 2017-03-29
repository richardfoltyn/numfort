
#include "numfort.h"

module numfort_stats_lm

    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_common_alloc
    use numfort_stats_core, only: normalize

    implicit none
    private

    interface ols
        module procedure ols_1d_real32, ols_1d_real64, ols_2d_real32, ols_2d_real64
    end interface

    interface ols_check_input
        module procedure ols_check_input_real32, ols_check_input_real64
    end interface

    interface pca
        module procedure pca_real32, pca_real64
    end interface

    interface pca_check_input
        module procedure pca_check_input_real32, pca_check_input_real64
    end interface

    interface pcr
        module procedure pcr_1d_real32, pcr_1d_real64, pcr_2d_real32, &
            pcr_2d_real64
    end interface

    interface pcr_check_input
        module procedure pcr_check_input_real32, pcr_check_input_real64
    end interface

    public :: ols, pca, pcr
contains


#define __PREC real64
#include "stats_lm_impl.F90"
#undef __PREC

#define __PREC real32
#include "stats_lm_impl.F90"
#undef __PREC

end module
