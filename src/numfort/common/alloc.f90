
module numfort_common_alloc
    !*  Module with helper routines to help management of dynamically allocated
    !   objects.

    use, intrinsic :: iso_fortran_env
    implicit none

    private

    public :: assert_alloc_ptr, assert_dealloc_ptr

    interface assert_alloc_ptr
        module procedure assert_alloc_ptr_1d_real32, assert_alloc_ptr_1d_real64
    end interface

    interface assert_dealloc_ptr
        module procedure assert_dealloc_ptr_1d_real32, assert_dealloc_ptr_1d_real64
    end interface

contains

! ------------------------------------------------------------------------------
! ASSERT_ALLOC_PTR routines

subroutine assert_alloc_ptr_1d_real32 (x, n, ptr_x)
    !*  ASSERT_ALLOC_PTR ensures that ptr_x points to allocated memory:
    !   If argument x is present, ptr_x points to x and no additional
    !   memory is allocated. If x is not present, then ptr_x points to a newly
    !   allocated block of memory of size n.

    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:), target, optional :: x
    integer, intent(in) :: n
    real (PREC), intent(out), dimension(:), pointer :: ptr_x

    if (present(x)) then
        if (size(x) >= n) then
            ptr_x => x
        else
            allocate (ptr_x(n))
        end if
    else
        allocate (ptr_x(n))
    end if
end subroutine

subroutine assert_alloc_ptr_1d_real64 (x, n, ptr_x)
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:), target, optional :: x
    integer, intent(in) :: n
    real (PREC), intent(out), dimension(:), pointer :: ptr_x

    if (present(x)) then
        if (size(x) >= n) then
            ptr_x => x
        else
            allocate (ptr_x(n))
        end if
    else
        allocate (ptr_x(n))
    end if
end subroutine

! ------------------------------------------------------------------------------
! ASSERT_DEALLOC_PTR routines

pure subroutine assert_dealloc_ptr_1d_real32 (x, ptr_x)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_x
    !   needs to be released only if argument x is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_x.

    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:), target, optional :: x
    real (PREC), intent(in out), dimension(:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine

pure subroutine assert_dealloc_ptr_1d_real64 (x, ptr_x)
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:), target, optional :: x
    real (PREC), intent(in out), dimension(:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine

end module
