

module numfort_stats_glmnet_real64

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: arange, split_uniform, argsort
    use numfort_common
    use numfort_common_alloc
    use numfort_common_input_checks
    use numfort_common_status
    use numfort_core, only: signum
    use numfort_linalg, only: gram

    use numfort_arrays, only: logspace, pack_indexed, unpack_indexed
    use numfort_stats_core
    use numfort_stats_data_helpers
    use numfort_stats_lm, only: predict

    use blas_interfaces, only: BLAS_AXPY => AXPY, BLAS_GEMV => GEMV, &
        BLAS_ASUM => ASUM, BLAS_NRM2 => NRM2, BLAS_IAMAX => IAMAX, &
        BLAS_GEMM => GEMM, BLAS_COPY => COPY, BLAS_GER => GER, BLAS_DOT => DOT

    use lapack_interfaces, only: LAPACK_GESDD => GESDD

    implicit none
    private


    public :: enet_config
    public :: enet_result

    public :: enet_path
    public :: enet_cv
    public :: enet_fit
    public :: enet_predict
    public :: enet_post_estim

    public :: predict
    public :: post_estim

    public :: create_alpha_grid
    public :: create_alpha_grid_cv

    public :: ridge

    integer, parameter :: PREC = real64

#include "include/glmnet_spec.F90"

    contains


#include "include/glmnet_impl.F90"


end module
