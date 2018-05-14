
#include <numfort.h>

module numfort_stats_lm

    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_common_alloc
    use numfort_stats_core, only: normalize, mean, std
    
    use blas_interfaces, only: BLAS_DOT => DOT, BLAS_GEMV => GEMV, BLAS_GEMM => GEMM
    use lapack_interfaces, only: LAPACK_GESVD => GESVD, LAPACK_GELSD => GELSD, &
        LAPACK_GESDD => GESDD

    implicit none
    private
    
    integer, public, parameter :: NF_STATS_LM_OLS = 1
    integer, public, parameter :: NF_STATS_LM_PCR = 2

    public :: ols, pca, pcr
    public :: finalize
    public :: post_estim

#include <numfort_real64.h>
#include "lm_data_spec.F90"
#include "stats_lm_spec.F90"

#include <numfort_real32.h>
#include "lm_data_spec.F90"
#include "stats_lm_spec.F90"

    contains

#include <numfort_real32.h>
#include "lm_data_impl.F90"
#include "stats_lm_impl.F90"

#include <numfort_real64.h>
#include "lm_data_impl.F90"
#include "stats_lm_impl.F90"

end module
