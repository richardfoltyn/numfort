
#include <numfort.h>

module numfort_stats_markov

    use, intrinsic :: iso_fortran_env

    use numfort_arrays_create, only: linspace
    use numfort_common_status
    use numfort_common_alloc
    use numfort_stats_dnorm, only: cdf, dnorm_real64
    use numfort_linalg, only: inv

    implicit none

    private
    public :: rouwenhorst
    public :: tauchen
    public :: ergodic_dist
    public :: moments

    interface rouwenhorst
        procedure rouwenhorst_real64
    end interface

    interface rouwenhorst_pad_matrix
        procedure rouwenhorst_pad_matrix_real64
    end interface

    interface tauchen
        procedure tauchen_real64
    end interface

    interface ergodic_dist
        procedure ergodic_dist_real64
    end interface

    interface moments
        procedure moments_real64
    end interface

    interface simulate
        procedure simulate_real64_int8, simulate_real64_int32
    end interface

    contains

#include <numfort_real64.h>
#include "markov_impl.F90"

#include <numfort_int8.h>
#include "markov_simulate_impl.F90"

#include <numfort_int32.h>
#include "markov_simulate_impl.F90"

end module
