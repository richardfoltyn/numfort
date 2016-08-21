#ifdef __INTEL_COMPILER
#define __SUPPORTS_PDT
#endif

#define __DEFAULT_REAL_KIND real64

#ifdef __SUPPORTS_PDT
#define __KIND_PARAM_DECL (PREC)
#define __REAL_KIND_DECL integer, kind :: PREC = __DEFAULT_REAL_KIND
#define __RWRK_KIND PREC
#else
#define __KIND_PARAM_DECL
#define __REAL_KIND_DECL
#define __RWRK_KIND __DEFAULT_REAL_KIND
#endif

module numfort_common

    use, intrinsic :: iso_fortran_env, only : real32, real64, int32

    implicit none
    private
    public :: workspace

    type, abstract :: workspace_base
        integer, dimension(:), allocatable :: iwrk
        logical, dimension(:), allocatable :: lwrk
        character (len=:), allocatable :: cwrk
    contains
        procedure (iface_assert_allocated), pass, deferred :: assert_allocated
    end type

    type, extends(workspace_base) :: workspace __KIND_PARAM_DECL
        __REAL_KIND_DECL
        real (__RWRK_KIND), dimension(:), allocatable :: rwrk
    contains
        procedure, pass, public :: assert_allocated => workspace_assert_allocated
    end type

    interface
        pure subroutine iface_assert_allocated (self, nrwrk, niwrk, ncwrk, nlwrk)
            import workspace_base
            class (workspace_base), intent(in out) :: self
            integer, intent(in), optional :: nrwrk, niwrk, ncwrk, nlwrk
        end subroutine
    end interface

    interface assert_allocated
        module procedure assert_allocated_real64, assert_allocated_real32, &
            assert_allocated_int, &
            assert_allocated_char, assert_allocated_bool
    end interface

contains

#ifdef __SUPPORTS_PDT
#define __WORKSPACE_TYPE workspace(PREC=__DEFAULT_REAL_KIND)
#else
#define __WORKSPACE_TYPE workspace
#endif

pure subroutine workspace_assert_allocated (self, nrwrk, niwrk, ncwrk, nlwrk)
    class (__WORKSPACE_TYPE), intent(in out) :: self
    integer, intent(in), optional :: nrwrk, niwrk, ncwrk, nlwrk

    call assert_allocated (self%rwrk, nrwrk)
    call assert_allocated (self%iwrk, niwrk)
    call assert_allocated (self%cwrk, ncwrk)
    call assert_allocated (self%lwrk, nlwrk)

end subroutine

pure subroutine assert_allocated_real64 (arr, n)
    real (real64), dimension(:) :: arr, tmp
    include "include/workspace_assert_allocated_impl.f90"
end subroutine

pure subroutine assert_allocated_real32 (arr, n)
    real (real32), dimension(:) :: arr, tmp
    include "include/workspace_assert_allocated_impl.f90"
end subroutine

pure subroutine assert_allocated_int(arr, n)
    integer, dimension(:) :: arr, tmp
    include "include/workspace_assert_allocated_impl.f90"
end subroutine

pure subroutine assert_allocated_char(arr, n)
    character (len=:), intent(in out), allocatable :: arr
    integer, intent(in), optional :: n
    character (len=:), allocatable :: tmp

    if (present(n)) then
        if (.not. allocated(arr)) then
            allocate (character (n) :: arr)
        else if (len(arr) < n) then
            allocate (character (n) :: tmp)
            call move_alloc (tmp, arr)
        end if
    end if
end subroutine

pure subroutine assert_allocated_bool(arr, n)
    logical, dimension(:) :: arr, tmp
    include "include/workspace_assert_allocated_impl.f90"
end subroutine

end module
