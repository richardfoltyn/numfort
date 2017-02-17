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

module numfort_common_workspace

    use, intrinsic :: iso_fortran_env, only : real32, real64, int32

    implicit none
    private
    public :: workspace

    integer, parameter :: SIZE_UNALLOCATED = -1

    type, abstract :: workspace_base
        integer, dimension(:), allocatable :: iwrk
        logical, dimension(:), allocatable :: lwrk
        character (len=:), allocatable :: cwrk
        integer :: nrwrk = SIZE_UNALLOCATED, niwrk = SIZE_UNALLOCATED, &
            ncwrk = SIZE_UNALLOCATED, nlwrk = SIZE_UNALLOCATED
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

    ! interface assert_allocated
    !     module procedure assert_allocated_real64, assert_allocated_real32, &
    !         assert_allocated_int, &
    !         assert_allocated_char, assert_allocated_bool
    ! end interface

contains

#ifdef __SUPPORTS_PDT
#define __WORKSPACE_TYPE workspace(PREC=__DEFAULT_REAL_KIND)
#else
#define __WORKSPACE_TYPE workspace
#endif

pure subroutine workspace_assert_allocated (self, nrwrk, niwrk, ncwrk, nlwrk)
    integer, parameter :: PREC = __DEFAULT_REAL_KIND
    class (__WORKSPACE_TYPE), intent(in out) :: self
    integer, intent(in), optional :: nrwrk, niwrk, ncwrk, nlwrk

    real (PREC), dimension(:), allocatable :: rtmp
    integer, dimension(:), allocatable :: itmp
    logical, dimension(:), allocatable :: ltmp
    character (:), allocatable :: ctmp

    ! Real working array
    if (present(nrwrk)) then
        if (nrwrk > 0) then
            if (.not. allocated(self%rwrk)) then
                allocate (self%rwrk(nrwrk))
            else if (size(self%rwrk) < nrwrk) then
                allocate (rtmp(nrwrk))
                rtmp(1:size(self%rwrk)) = self%rwrk
                call move_alloc (rtmp, self%rwrk)
            end if
            self%nrwrk = nrwrk
        end if
    end if

    ! Integer working array
    if (present(niwrk)) then
        if (niwrk > 0) then
            if (.not. allocated(self%iwrk)) then
                allocate (self%iwrk(niwrk))
            else if (size(self%iwrk) < niwrk) then
                allocate (itmp(niwrk))
                itmp(1:size(self%iwrk)) = self%iwrk
                call move_alloc (itmp, self%iwrk)
            end if
            self%niwrk = niwrk
        end if
    end if

    ! Logical working array
    if (present(nlwrk)) then
        if (nlwrk > 0) then
            if (.not. allocated(self%lwrk)) then
                allocate (self%lwrk(nlwrk))
            else if (size(self%lwrk) < nlwrk) then
                allocate (ltmp(nlwrk))
                ltmp(1:size(self%lwrk)) = self%lwrk
                call move_alloc (ltmp, self%lwrk)
            end if
            self%nlwrk = nlwrk
        end if
    end if

    if (present(ncwrk)) then
        if (ncwrk > 0) then
            if (.not. allocated(self%cwrk)) then
                allocate (character (ncwrk) :: self%cwrk)
            else if (len(self%cwrk) < ncwrk) then
                allocate (character (ncwrk) :: ctmp)
                call move_alloc (ctmp, self%cwrk)
            end if
            self%ncwrk = ncwrk
        end if
    end if

    ! call assert_allocated (self%rwrk, self%nrwrk, nrwrk)
    ! call assert_allocated (self%iwrk, self%niwrk, niwrk)
    ! call assert_allocated (self%cwrk, self%ncwrk, ncwrk)
    ! call assert_allocated (self%lwrk, self%nlwrk, nlwrk)

end subroutine
!
! pure subroutine assert_allocated_real64 (arr, n_attr, n)
!     real (real64), dimension(:) :: arr, tmp
!     include "include/workspace_assert_allocated_impl.f90"
! end subroutine
!
! pure subroutine assert_allocated_real32 (arr, n_attr, n)
!     real (real32), dimension(:) :: arr, tmp
!     include "include/workspace_assert_allocated_impl.f90"
! end subroutine
!
! pure subroutine assert_allocated_int(arr, n_attr, n)
!     integer, dimension(:) :: arr, tmp
!     include "include/workspace_assert_allocated_impl.f90"
! end subroutine
!
! pure subroutine assert_allocated_char(arr, n_attr, n)
!     character (len=:), intent(in out), allocatable :: arr
!     integer, intent(in out) :: n_attr
!     integer, intent(in), optional :: n
!     character (len=:), allocatable :: tmp
!
!     if (present(n)) then
!         if (n > 0) then
!             if (.not. allocated(arr)) then
!                 allocate (character (n) :: arr)
!                 n_attr = n
!             else if (len(arr) < n) then
!                 allocate (character (n) :: tmp)
!                 call move_alloc (tmp, arr)
!                 n_attr = n
!             end if
!         end if
!     end if
! end subroutine
!
! pure subroutine assert_allocated_bool(arr, n_attr, n)
!     logical, dimension(:) :: arr, tmp
!     include "include/workspace_assert_allocated_impl.f90"
! end subroutine

end module
