
#include <numfort.h>


module numfort_linalg_core
    !*  Module contains some core linear algebra routine such as
    !   computing matrix inverses and determinants.

    use, intrinsic :: iso_fortran_env
    use numfort_common_status
    use numfort_common, only: shape_equal
    use numfort_common_workspace
    use lapack_interfaces, only: LAPACK_GETRF => GETRF, LAPACK_GETRI => GETRI

    implicit none

    private
    public :: inv
    public :: inv_work_query
    public :: det

#include <numfort_real32.h>
#include "linalg_core_spec.F90"

#include <numfort_real64.h>
#include "linalg_core_spec.F90"

    contains


#include <numfort_real32.h>
#include "linalg_core_impl.F90"

#include <numfort_real64.h>
#include "linalg_core_impl.F90"

end module
