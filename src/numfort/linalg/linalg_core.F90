
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

    interface inv_work_query
        procedure inv_work_query_dims
    end interface

    contains


#include <numfort_real32.h>
#include "linalg_core_impl.F90"

#include <numfort_real64.h>
#include "linalg_core_impl.F90"

! ------------------------------------------------------------------------------
! Precision-independent routines

subroutine inv_work_query_dims (n, lwork, liwork, status)
    !*  INV_WORK_QUERY returns the optimal workspace sizes for real and
    !   integer working arrays needed by the INV routine.
    integer, parameter :: PREC = __PREC
    integer, intent(in) :: n
        !*  Number of rows (columns) of square matrix to be inverted.
    integer, intent(out) :: lwork
        !*  Optimal REAL working array size
    integer, intent(out) :: liwork
        !*  INTEGER working array size
    type (status_t), intent(out), optional :: status
        !*  Exit status.

    type (status_t) :: lstatus
    integer :: lda, info
    integer, dimension(1) :: idummy1d
    real (PREC), dimension(1) :: rdummy1d
    real (PREC), dimension(1,1) :: rdummy2d

    lstatus = NF_STATUS_OK

    liwork = -1

    ! Workspace query for GETRI
    lda = n
    lwork = -1
    call LAPACK_GETRI (n, rdummy2d, lda, idummy1d, rdummy1d, lwork, info)
    if (info /= 0) then
        lstatus = NF_STATUS_INVALID_STATE
        goto 100
    end if

    ! Optimal array size stored in first element of RWORK1
    lwork = int(rdummy1d(1))
    liwork = n

100 continue

    if (present(status)) status = lstatus

end subroutine


end module
