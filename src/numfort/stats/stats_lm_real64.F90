

module numfort_stats_lm_real64

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_alloc
    use numfort_common_input_checks
    use numfort_stats_core
    use numfort_stats_lm_common
    use numfort_stats_data_helpers
    
    use blas_interfaces, only: BLAS_DOT => DOT, BLAS_GEMV => GEMV, &
        BLAS_GEMM => GEMM, BLAS_NRM2 => NRM2, BLAS_COPY => COPY, BLAS_AXPY => AXPY
    use lapack_interfaces, only: LAPACK_GESVD => GESVD, LAPACK_GELSD => GELSD, &
        LAPACK_GESDD => GESDD

    implicit none
    private

    public :: lm_result
    public :: lm_config

    public :: lm_result_reset

    public :: ols
    public :: pca
    public :: pcr
    public :: predict

    public :: post_estim


    integer, parameter :: PREC = real64


#include "stats_lm_spec.F90"

    contains

#include "stats_lm_impl.F90"

end module
