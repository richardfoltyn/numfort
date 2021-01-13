

module numfort_stats_np_real64
    !*  Module implementing some core non-parametric functionality
    !   for double-precision data

    use, intrinsic :: iso_fortran_env

    use numfort_core, only: PI => PI_real64

    implicit none

    private

    public :: gaussian_kde

    interface gaussian_kde
        procedure gaussian_kde
    end interface

    integer, parameter :: PREC = real64


    contains


#include "stats_np_impl.F90"



end module
