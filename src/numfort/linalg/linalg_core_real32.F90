

module numfort_linalg_core_real32
    !*  Module contains some core linear algebra routine such as
    !   computing matrix inverses and determinants.

    use, intrinsic :: iso_fortran_env

    use numfort_common_status
    use numfort_common, only: shape_equal
    use numfort_common_workspace, workspace => workspace_real32
    use lapack_interfaces, only: LAPACK_GETRF => GETRF, LAPACK_GETRI => GETRI

    implicit none

    private

    public :: inv
    public :: inv_work_query
    public :: det

    interface inv
        procedure inv
    end interface

    interface inv_work_query
        procedure inv_work_query_dims, inv_work_query_matrix
    end interface

    interface det
        procedure det
    end interface

    integer, parameter :: PREC = real32

    contains

#include "include/linalg_core_impl.F90"


end module
