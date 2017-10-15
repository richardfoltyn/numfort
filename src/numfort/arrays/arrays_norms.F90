
#include "numfort.h"

module numfort_arrays_norms

    use, intrinsic :: iso_fortran_env
    use numfort_arrays_shape, only: shape_equal
    implicit none
    private

    public :: norm_sup_diff

    interface norm_sup_diff
        module procedure norm_sup_diff_3d_real64
    end interface

    interface norm_sup_diff
        module procedure norm_sup_diff_4d_real64
    end interface
contains


#define __PREC real64
#include "norms_impl.F90"
#undef __PREC


end module
