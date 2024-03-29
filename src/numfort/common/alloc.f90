
module numfort_common_alloc
    !*  Module with helper routines to help management of dynamically allocated
    !   objects.

    use, intrinsic :: iso_fortran_env
    implicit none

    private

    public :: assert_alloc_ptr, assert_dealloc_ptr

    interface assert_alloc_ptr
        procedure assert_alloc_ptr_1d_real32, assert_alloc_ptr_1d_real64, &
            assert_alloc_ptr_1d_int32, assert_alloc_ptr_1d_int64
    end interface

    interface assert_dealloc_ptr
        procedure assert_dealloc_ptr_1d_logical, &
            assert_dealloc_ptr_1d_int32, assert_dealloc_ptr_1d_int64, &
            assert_dealloc_ptr_1d_real32, assert_dealloc_ptr_1d_real64
    end interface

    interface assert_dealloc_ptr
        procedure assert_dealloc_ptr_2d_logical, &
            assert_dealloc_ptr_2d_real32, assert_dealloc_ptr_2d_real64
    end interface

contains

! ------------------------------------------------------------------------------
! ASSERT_ALLOC_PTR routines

pure subroutine assert_alloc_ptr_1d_real32 (x, n, ptr_x)
    !*  ASSERT_ALLOC_PTR ensures that ptr_x points to allocated memory:
    !   If argument x is present, ptr_x points to x and no additional
    !   memory is allocated. If x is not present or of size zero, then ptr_x
    !   points to a newly allocated block of memory of size n.
    integer, parameter :: PREC = real32

    real (PREC), intent(inout), dimension(:), target, optional :: x
        !*  Target array
    integer, intent(in) :: n
        !*  Minimal array size
    real (PREC), intent(out), dimension(:), pointer :: ptr_x

    nullify (ptr_x)

    if (present(x)) then
        if (size(x) == 0) then
            ! Need to handle size-zero arrays separately, as
            ! associated (ptr_x, x) will return .FALSE. if ptr_x => x
            ! but x is of size zero.
            ! Then assert_dealloc_ptr will attempt to free memory that was
            ! never allocated!
            allocate (ptr_x(0))
        else if (size(x) >= n) then
            ptr_x => x
        else
            allocate (ptr_x(n))
        end if
    else
        allocate (ptr_x(n))
    end if
end subroutine



pure subroutine assert_alloc_ptr_1d_real64 (x, n, ptr_x)
    integer, parameter :: PREC = real64

    real (PREC), intent(inout), dimension(:), target, optional :: x
    integer, intent(in) :: n
    real (PREC), intent(out), dimension(:), pointer :: ptr_x

    nullify (ptr_x)

    if (present(x)) then
        if (size(x) == 0) then
            ! Need to handle size-zero arrays separately, as
            ! associated (ptr_x, x) will return .FALSE. if ptr_x => x
            ! but x is of size zero.
            ! Then assert_dealloc_ptr will attempt to free memory that was
            ! never allocated!
            allocate (ptr_x(0))
        else if (size(x) >= n) then
            ptr_x => x
        else
            allocate (ptr_x(n))
        end if
    else
        allocate (ptr_x(n))
    end if
end subroutine



pure subroutine assert_alloc_ptr_1d_int32 (x, n, ptr_x)
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(inout), dimension(:), target, optional :: x
    integer, intent(in) :: n
    integer (INTSIZE), intent(out), dimension(:), pointer :: ptr_x

    nullify (ptr_x)

    if (present(x)) then
        if (size(x) == 0) then
            ! Need to handle size-zero arrays separately, as
            ! associated (ptr_x, x) will return .FALSE. if ptr_x => x
            ! but x is of size zero.
            ! Then assert_dealloc_ptr will attempt to free memory that was
            ! never allocated!
            allocate (ptr_x(0))
        else if (size(x) >= n) then
            ptr_x => x
        else
            allocate (ptr_x(n))
        end if
    else
        allocate (ptr_x(n))
    end if
end subroutine



pure subroutine assert_alloc_ptr_1d_int64 (x, n, ptr_x)
    integer, parameter :: INTSIZE = int64

    integer (INTSIZE), intent(inout), dimension(:), target, optional :: x
    integer, intent(in) :: n
    integer (INTSIZE), intent(out), dimension(:), pointer :: ptr_x

    nullify (ptr_x)

    if (present(x)) then
        if (size(x) == 0) then
            ! Need to handle size-zero arrays separately, as
            ! associated (ptr_x, x) will return .FALSE. if ptr_x => x
            ! but x is of size zero.
            ! Then assert_dealloc_ptr will attempt to free memory that was
            ! never allocated!
            allocate (ptr_x(0))
        else if (size(x) >= n) then
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

pure subroutine assert_dealloc_ptr_1d_logical (x, ptr_x)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_x
    !   needs to be released only if argument x is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_x.

    logical, intent(in), dimension(:), target, optional :: x
    logical, intent(inout), dimension(:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine


pure subroutine assert_dealloc_ptr_1d_int32 (x, ptr_x)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_x
    !   needs to be released only if argument x is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_x.

    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:), target, optional :: x
    integer (INTSIZE), intent(inout), dimension(:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine

pure subroutine assert_dealloc_ptr_1d_int64 (x, ptr_x)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_x
    !   needs to be released only if argument x is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_x.

    integer, parameter :: INTSIZE = int64

    integer (INTSIZE), intent(in), dimension(:), target, optional :: x
    integer (INTSIZE), intent(inout), dimension(:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine

pure subroutine assert_dealloc_ptr_1d_real32 (x, ptr_x)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_x
    !   needs to be released only if argument x is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_x.

    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:), target, optional :: x
    real (PREC), intent(inout), dimension(:), pointer :: ptr_x

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
    real (PREC), intent(inout), dimension(:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine

pure subroutine assert_dealloc_ptr_2d_logical (x, ptr_x)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_x
    !   needs to be released only if argument x is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_x.

    logical, intent(in), dimension(:,:), target, optional :: x
    logical, intent(inout), dimension(:,:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine

pure subroutine assert_dealloc_ptr_2d_real32 (x, ptr_x)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_x
    !   needs to be released only if argument x is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_x.

    integer, parameter :: PREC = real32

    real (PREC), intent(in), dimension(:,:), target, optional :: x
    real (PREC), intent(inout), dimension(:,:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine

pure subroutine assert_dealloc_ptr_2d_real64 (x, ptr_x)
    integer, parameter :: PREC = real64

    real (PREC), intent(in), dimension(:,:), target, optional :: x
    real (PREC), intent(inout), dimension(:,:), pointer :: ptr_x

    if (associated(ptr_x)) then
        if (.not. present(x)) then
            deallocate (ptr_x)
        else if (.not. associated(ptr_x, x)) then
            deallocate (ptr_x)
        end if
    end if
end subroutine


end module
