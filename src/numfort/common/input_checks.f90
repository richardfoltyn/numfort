module numfort_common_input_checks
    !*  Module for helper functions that perform input validation.

    use, intrinsic :: iso_fortran_env
    use numfort_common_kinds
    use numfort_common_status
    implicit none

    private

    public :: check_positive, check_nonneg, check_enum

    interface check_positive
        module procedure check_positive_real64, check_positive_int32
    end interface

    interface check_nonneg
        module procedure check_nonneg_real64, check_nonneg_int32
    end interface

contains

! ------------------------------------------------------------------------------
! CHECK_POSITIVE routines

pure subroutine check_positive_real64 (dummy, val, name, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val

    include "include/check_positive.f90"
end subroutine

pure subroutine check_positive_int32 (dummy, val, name, status, msg)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: dummy
    integer (INTSIZE), intent(in), optional :: val

    include "include/check_positive.f90"
end subroutine

! ------------------------------------------------------------------------------
! CHECK_NONNEG implementations

pure subroutine check_nonneg_real64 (dummy, val, name, status, msg)
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: dummy
    real (PREC), intent(in), optional :: val

    include "include/check_nonneg.f90"
end subroutine

pure subroutine check_nonneg_int32 (dummy, val, name, status, msg)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: dummy
    integer (INTSIZE), intent(in), optional :: val

    include "include/check_nonneg.f90"
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
    type (status_t), intent(in out) :: status
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
    type (status_t), intent(in out) :: status
    character (*), intent(out), optional :: msg

    status = NF_STATUS_OK
    if (present(msg)) msg = ""
end subroutine

end module
