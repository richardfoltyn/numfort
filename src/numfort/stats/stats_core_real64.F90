

module numfort_stats_core_real64

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: argsort
    use numfort_common
    use numfort_interpolate
    use numfort_stats_core_common

    use blas_interfaces, only: BLAS_GER => GER, BLAS_GEMM => GEMM

    implicit none

    integer, parameter :: PREC = real64

    private

    public :: mean
    public :: std
    public :: cov
    public :: corrcoef
    public :: normalize
    public :: quantile

    interface mean
        procedure mean_1d, mean_2d
    end interface

    interface std
        procedure std_1d, std_2d
    end interface

    interface normalize
        procedure normalize_2d
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
