

module numfort_stats_markov_simulate_real64_int8

    use, intrinsic :: iso_fortran_env

    use numfort_common_status
    use numfort_stats_markov_real64

    implicit none
    private

    public :: markov_simulate
    public :: markov_simulate_advanced

    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int8

    interface markov_simulate
        procedure simulate
    end interface

    interface markov_simulate_advanced
        procedure simulate_advanced
    end interface


    contains

#include "markov_simulate_impl.F90"


end module
