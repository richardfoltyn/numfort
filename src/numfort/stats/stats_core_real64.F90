

module numfort_stats_core_real64

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: argsort
    use numfort_common_alloc
    use numfort_common_input_checks
    use numfort_common_status
    use numfort_interpolate
    use numfort_stats_core_common

    use blas_interfaces, only: BLAS_GER => GER, BLAS_GEMM => GEMM

    implicit none

    integer, parameter :: PREC = real64

    private

    public :: mean
    public :: mean_impl
    public :: std
    public :: std_impl
    public :: normalize
    public :: standardize
    public :: standardize_impl
    public :: cov
    public :: corrcoef
    public :: quantile

    interface mean
        procedure mean_1d, mean_2d
    end interface

    interface mean_impl
        procedure mean_impl_1d, mean_impl_2d
    end interface

    interface std
        procedure std_1d, std_2d
    end interface

    interface std_impl
        procedure std_impl_1d, std_impl_2d
    end interface

    interface normalize
        procedure normalize_2d
    end interface

    interface standardize
        procedure standardize_1d, standardize_2d
    end interface

    interface standardize_impl
        procedure standardize_impl_1d, standardize_impl_2d
    end interface

    interface quantile
        procedure quantile_dispatch, quantile_dispatch_scalar
    end interface

    interface cov
        procedure cov
    end interface

    interface corrcoef
        procedure corrcoef
    end interface

    contains

#include "stats_core_impl.F90"

end module
