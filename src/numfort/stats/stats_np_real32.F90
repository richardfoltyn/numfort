

module numfort_stats_np_real32
    !*  Module implementing some core non-parametric functionality
    !   for double-precision data

    use, intrinsic :: iso_fortran_env

    use numfort_core, only: PI => PI_real32

    implicit none

    private

    public :: gaussian_kde

    interface gaussian_kde
        procedure gaussian_kde, gaussian_kde_1d
    end interface

    integer, parameter :: PREC = real32


    contains


#include "stats_np_impl.F90"



end module
