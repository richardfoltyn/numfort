

module numfort_stats_ineq_real32
    !*  Module implements inequality measures such as the Gini coefficient

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays, only: argsort

    implicit none
    private

    public :: gini

    integer, parameter :: PREC = real32

    interface gini
        procedure gini_impl
    end interface

    contains

#include "stats_ineq_impl.F90"

end module