

module numfort_stats_data_helpers_real64

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays, only: arange, unique
    use numfort_common_alloc
    use numfort_common_cond_alloc
    use numfort_common_status
    use numfort_common_input_checks
    use numfort_stats_core

    use blas_interfaces, only: BLAS_GEMV => GEMV, BLAS_GEMM => GEMM

    implicit none
    private

    public :: transform_regressors
    public :: replay_transform
    public :: has_all_vars
    public :: get_regressor_dims
    public :: extract_block
    public :: extract_block_alloc
    public :: random_sample

    integer, parameter :: PREC = real64

    interface transform_regressors
        procedure transform_regr, transform_regr_in_place
    end interface

    interface replay_transform
        procedure replay_transform, replay_transform_in_place
    end interface

    interface get_regressor_dims
        procedure get_regressor_dims
    end interface

    interface has_all_vars
        procedure has_all_vars
    end interface

    interface extract_block
        procedure extract_block_1d, extract_block_2d
    end interface

    interface extract_block_alloc
        procedure extract_block_alloc_1d, extract_block_alloc_2d
    end interface

    interface random_sample
        procedure random_sample_1d, random_sample_2d
    end interface

    contains


#include "include/data_helpers_impl.F90"

end module
