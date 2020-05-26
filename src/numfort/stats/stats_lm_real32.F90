

module numfort_stats_lm_real32

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_alloc
    use numfort_common_input_checks
    use numfort_stats_core_real32, only: normalize, mean, std
    use numfort_stats_lm_common

    use blas_interfaces, only: BLAS_DOT => DOT, BLAS_GEMV => GEMV, BLAS_GEMM => GEMM
    use lapack_interfaces, only: LAPACK_GESVD => GESVD, LAPACK_GELSD => GELSD, &
        LAPACK_GESDD => GESDD

    implicit none
    private

    public :: ols, pca, pcr
    public :: finalize
    public :: ridge
    public :: post_estim

    public :: lm_data

    integer, parameter :: PREC = real32

    interface ols
        procedure ols_1d, ols_2d
    end interface

    interface pca
        procedure pca
    end interface

    interface pcr
        procedure pcr_1d, pcr_2d
    end interface

    interface pcr
        procedure pcr_pca_1d, pcr_pca_2d
    end interface

    interface ridge
        procedure ridge_1d, ridge_2d
    end interface

    interface post_estim
        procedure lm_post_estim
    end interface

#include "lm_data_spec.F90"

    contains

#include "lm_data_impl.F90"
#include "stats_lm_impl.F90"

end module
