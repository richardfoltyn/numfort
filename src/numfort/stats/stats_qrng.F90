

#include "numfort.h"

module numfort_stats_qrng
    !*  Module implements generation of quasi-random sequences (such as
    !   Sobol sequences).

    use, intrinsic :: iso_fortran_env
    use sobol_direction_num
    use numfort_common

    implicit none

    private

    public :: sobol_state
    public :: sobol_init, sobol_reset, sobol_next

    ! Do not use sign bit so all numbers remain positive!
    integer, parameter :: SOBOL_MAX_BITS = 63

    type :: sobol_state
        private
        integer (int64) :: ndim = 1
            !*  Dimension of Sobol sequence
        integer (int64) :: idx = 0
            !*  Position of last sampled point in Sobol sequence
        integer (int64), dimension(:), allocatable :: x
            !*  Element of Sobol sequence at position idx.
        integer (int64), dimension(:,:), allocatable :: v
            !*  Direction numbers, scaled by 2 ** SOBOL_MAX_BITS; each
            !   column contains up to SOBOL_MAX_BITS direction numbers for
            !   the corresponding dimension.
    end type

    interface sobol_init_check_input
        module procedure sobol_init_check_input_int32, &
            sobol_init_check_input_int64
    end interface

    interface sobol_init
        module procedure sobol_init_scalar_int32, sobol_init_array_int32, &
            sobol_init_scalar_int64, sobol_init_array_int64
    end interface

    interface sobol_next
        module procedure sobol_next_scalar_real32, sobol_next_array_real32, &
            sobol_next_scalar_real64, sobol_next_array_real64
    end interface

contains

pure subroutine sobol_reset (self)
    type (sobol_state), intent(in out) :: self

    self%ndim = 1
    self%idx = 0
    if (allocated(self%x)) deallocate(self%x)
    if (allocated(self%v)) deallocate(self%v)
end subroutine

#define __INTSIZE int32
#include "sobol_init_impl.F90"
#undef __INTSIZE

#define __INTSIZE int64
#include "sobol_init_impl.F90"
#undef __INTSIZE

#define __PREC real32
#include "sobol_eval_impl.F90"
#undef __PREC

#define __PREC real64
#include "sobol_eval_impl.F90"
#undef __PREC



end module
