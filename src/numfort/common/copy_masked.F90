

module numfort_common_copy_masked

    use, intrinsic :: iso_fortran_env

    use numfort_common_status
    use numfort_common_shape
    use numfort_common_cond_alloc

    implicit none

    private

    public :: copy

    interface copy
        procedure copy_masked_1d_real32, copy_masked_1d_real64
    end interface

    interface copy
        procedure copy_masked_2d_real32, copy_masked_2d_real64
    end interface

    contains


pure subroutine copy_masked_1d_real32 (src, dst, mask, status)
    !*  COPY_MASKED_1D is a wrapper around the intrinsic PACK routine
    !   which does not create temporary arrays when some of the objects
    !   involved have a POINTER or TARGET attribute.
    integer, parameter :: PREC = real32

#include "copy_masked_1d_impl.F90"
end subroutine

pure subroutine copy_masked_1d_real64 (src, dst, mask, status)
    !*  COPY_MASKED_1D is a wrapper around the intrinsic PACK routine
    !   which does not create temporary arrays when some of the objects
    !   involved have a POINTER or TARGET attribute.
    integer, parameter :: PREC = real64

#include "copy_masked_1d_impl.F90"
end subroutine


pure subroutine copy_masked_2d_real32 (src, dst, mask, dim, status)
    !*  COPY_MASKED_2D copies a subject of rows or columns in SRC into DST
    !   where values to be copied are identified by MASK.
    integer, parameter :: PREC = real32

#include "copy_masked_2d_impl.F90"
end subroutine

pure subroutine copy_masked_2d_real64 (src, dst, mask, dim, status)
    integer, parameter :: PREC = real64
    !*  COPY_MASKED_2D copies a subject of rows or columns in SRC into DST
    !   where values to be copied are identified by MASK.

#include "copy_masked_2d_impl.F90"
end subroutine



end module
