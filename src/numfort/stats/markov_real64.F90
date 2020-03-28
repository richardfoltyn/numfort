

module numfort_stats_markov_real64

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: arange, linspace, setdiff
    use numfort_common
    use numfort_common_alloc
    use numfort_common_input_checks
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_stats_dnorm, only: cdf, dnorm => dnorm_real64
    use numfort_linalg, only: inv, inv_work_query

    use blas_interfaces, only: BLAS_GEMV => GEMV, BLAS_COPY => COPY, &
        BLAS_AXPY => AXPY

    implicit none

    private

    public :: rouwenhorst
    public :: tauchen
    public :: markov_ergodic_dist
    public :: markov_moments
    public :: is_trans_matrix
    public :: truncate_trans_matrix

    integer, parameter :: PREC = real64

    interface rouwenhorst
        procedure rouwenhorst
    end interface

    interface tauchen
        procedure tauchen
    end interface

    interface markov_ergodic_dist
        procedure ergodic_dist
    end interface

    interface markov_moments
        procedure moments
    end interface

    interface is_trans_matrix
        procedure is_trans_matrix
    end interface

    interface truncate_trans_matrix
        procedure truncate_trans_matrix
    end interface

    contains

#include "markov_impl.F90"

end module
