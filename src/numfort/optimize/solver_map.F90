
#include <numfort.h>

module numfort_optimize_solver_map

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays
    use numfort_common_status
    use numfort_common_input_checks

    implicit none
    private


    public :: solver_map
    public :: solver_map_init
    public :: solver_map_eval
    public :: solver_map_eval_inverse

    integer, parameter :: TRANSFORM_LINEAR = 1
    integer, parameter :: TRANSFORM_EXP = 2
    integer, parameter :: TRANSFORM_NEG_EXP = 3
    integer, parameter :: TRANSFORM_LOGISTIC = 4

    type :: solver_map
        private
        integer :: transform = TRANSFORM_LINEAR
        real (real64) :: lb
        real (real64) :: ub
        real (real64) :: scale = 1.0
    end type


#include <numfort_real32.h>
#include "solver_map_spec.F90"

#include <numfort_real64.h>
#include "solver_map_spec.F90"


    contains


#include <numfort_real32.h>
#include "solver_map_impl.F90"

#include <numfort_real64.h>
#include "solver_map_impl.F90"

end module
