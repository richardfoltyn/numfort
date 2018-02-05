
#include <numfort.h>

module numfort_stats_markov

    use, intrinsic :: iso_fortran_env

    use numfort_arrays_create, only: linspace
    use numfort_common_status
    use numfort_common_alloc
    use numfort_common_swap
    use numfort_stats_dnorm, only: cdf, dnorm_real64
    use numfort_linalg, only: inv

    use blas_interfaces, only: GEMV, COPY, AXPY

    implicit none

    private
    public :: rouwenhorst
    public :: tauchen
    public :: markov_ergodic_dist
    public :: markov_moments
    public :: markov_simulate
    public :: is_trans_matrix

    interface rouwenhorst
        procedure rouwenhorst_real64
    end interface

    interface rouwenhorst_pad_matrix
        procedure rouwenhorst_pad_matrix_real64
    end interface

    interface tauchen
        procedure tauchen_real64
    end interface

    interface markov_ergodic_dist
        procedure ergodic_dist_real64
    end interface

    interface markov_moments
        procedure moments_real64
    end interface

    interface markov_simulate
        procedure simulate_real64_int8, simulate_real64_int32
    end interface

    interface is_trans_matrix
        procedure is_trans_matrix_real64
    end interface
    
    interface markov_approx_input_checks
        procedure markov_approx_input_checks_real64
    end interface

    contains

#include <numfort_real64.h>
#include "markov_impl.F90"

#include <numfort_int8.h>
#include "markov_simulate_impl.F90"

#include <numfort_int32.h>
#include "markov_simulate_impl.F90"

end module
