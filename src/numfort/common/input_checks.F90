module numfort_common_input_checks
    !*  Module for helper functions that perform input validation.

    use, intrinsic :: iso_fortran_env
    use numfort_common_kinds
    use numfort_common_status
    implicit none

    private

    public :: check_cond
    public :: check_positive
    public :: check_nonneg
    public :: check_nonzero
    public :: check_enum
    public :: check_range

    interface check_positive
        procedure check_positive_real32, check_positive_real64, &
            check_positive_int32
    end interface

    interface check_nonneg
        procedure check_nonneg_real32, check_nonneg_real64, check_nonneg_int32
    end interface

    interface check_nonzero
        procedure check_nonzero_real32, check_nonzero_real64
    end interface

    interface check_range
        procedure check_range_real32, check_range_real64
    end interface

contains


! ------------------------------------------------------------------------------
! Check generic logical condition

pure subroutine check_cond (cond, prefix, lmsg, status, msg)
    logical, intent(in) :: cond
    character (*), intent(in) :: prefix
    character (*), intent(in) :: lmsg
    type (status_t), intent(out) :: status
    character (*), intent(out), optional :: msg

    status = NF_STATUS_OK
    if (.not. cond) then
        status = NF_STATUS_INVALID_ARG
        if (present(msg)) then
            if (len_trim(prefix) > 0) then
                msg = prefix // ': ' // lmsg
            else
                msg = lmsg
            end if
        end if
    end if
end subroutine


! ------------------------------------------------------------------------------
! CHECK_POSITIVE routines

pure subroutine check_positive_real64 (dummy, val, name, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val

#include "check_positive.f90"
end subroutine

pure subroutine check_positive_real32 (dummy, val, name, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val

#include "check_positive.f90"
end subroutine

pure subroutine check_positive_int32 (dummy, val, name, status, msg)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: dummy
    integer (INTSIZE), intent(in), optional :: val

#include "check_positive.f90"
end subroutine

! ------------------------------------------------------------------------------
! CHECK_NONNEG implementations

pure subroutine check_nonneg_real32 (dummy, val, name, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val

#include "check_nonneg.f90"
end subroutine

pure subroutine check_nonneg_real64 (dummy, val, name, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val

#include "check_nonneg.f90"
end subroutine

pure subroutine check_nonneg_int32 (dummy, val, name, status, msg)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: dummy
    integer (INTSIZE), intent(in), optional :: val

#include "check_nonneg.f90"
end subroutine


! ------------------------------------------------------------------------------
! CHECK_NONZERO implementations

pure subroutine check_nonzero_real32 (dummy, val, name, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val

#include "check_nonzero.f90"
end subroutine

pure subroutine check_nonzero_real64 (dummy, val, name, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val

#include "check_nonzero.f90"
end subroutine


! ------------------------------------------------------------------------------
! CHECK_RANGE implementations

pure subroutine check_range_real32 (dummy, val, lb, ub, name, status, msg)
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val
    real (PREC), intent(in), optional :: lb, ub

#include "check_range.f90"
end subroutine

pure subroutine check_range_real64 (dummy, val, lb, ub, name, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val
    real (PREC), intent(in), optional :: lb, ub

#include "check_range.f90"
end subroutine

! ------------------------------------------------------------------------------

pure subroutine check_enum (val, valid, name, status, msg)
    !*  CHECK_ENUM verifies that the value of the optional argument, if present,
    !   is an element of the set of valid values.

    integer (NF_ENUM_KIND), intent(in), optional :: val
        !*  User-provided argument
    integer (NF_ENUM_KIND), intent(in), dimension(:) :: valid
        !*  Array containing the set of valid enum values.
    character (*), intent(in), optional :: name
        !*  Optional argument name, used in error message.
    type (status_t), intent(inout) :: status
    character (*), intent(out), optional :: msg

    integer :: i

    call clear_status (status, msg)

    if (present(val)) then
        do i = 1, size(valid)
            if (val == valid(i)) return
        end do

        ! At this point no matching value was found
        status = NF_STATUS_INVALID_ARG
        if (present(msg)) then
            if (present(name)) then
                msg = "Argument '" // trim(name) // "': invalid value"
            else
                msg = "Invalid value in argument"
            end if
        end if
   end if
end subroutine

! ------------------------------------------------------------------------------
! HELPER ROUTINES

pure subroutine clear_status (status, msg)
    type (status_t), intent(inout) :: status
    character (*), intent(out), optional :: msg

    status = NF_STATUS_OK
    if (present(msg)) msg = ""
end subroutine

end module
