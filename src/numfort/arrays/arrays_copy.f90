module numfort_arrays_copy

    use, intrinsic :: iso_fortran_env

    use numfort_common_shape

    implicit none

    private
    public :: copy_alloc


    interface copy_alloc
        module procedure copy_alloc_1d_real32, copy_alloc_1d_real64, &
            copy_alloc_1d_int32, copy_alloc_1d_int64
    end interface

    interface copy_alloc
        module procedure copy_alloc_2d_real32, copy_alloc_2d_real64, &
            copy_alloc_2d_int32, copy_alloc_2d_int64
    end interface

    interface copy_alloc
        module procedure copy_alloc_3d_real32, copy_alloc_3d_real64, &
            copy_alloc_3d_int32, copy_alloc_3d_int64
    end interface

    interface copy_alloc
        module procedure copy_alloc_4d_real32, copy_alloc_4d_real64, &
            copy_alloc_4d_int32, copy_alloc_4d_int64
    end interface

    contains


pure subroutine copy_alloc_1d_real64 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: PREC = real64
    integer, parameter :: ND = 1

    real (PREC), intent(in), dimension(:), optional :: src
    real (PREC), intent(out), dimension(:), allocatable :: dst
    include "include/copy_alloc_1d_impl.f90"
end subroutine

pure subroutine copy_alloc_1d_real32 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: PREC = real32
    integer, parameter :: ND = 1

    real (PREC), intent(in), dimension(:), optional :: src
    real (PREC), intent(out), dimension(:), allocatable :: dst
    include "include/copy_alloc_1d_impl.f90"
end subroutine

pure subroutine copy_alloc_1d_int32 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: INTSIZE = int32
    integer, parameter :: ND = 1

    integer (INTSIZE), intent(in), dimension(:), optional :: src
    integer (INTSIZE), intent(out), dimension(:), allocatable :: dst
    include "include/copy_alloc_1d_impl.f90"
end subroutine

pure subroutine copy_alloc_1d_int64 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: INTSIZE = int64
    integer, parameter :: ND = 1

    integer (INTSIZE), intent(in), dimension(:), optional :: src
    integer (INTSIZE), intent(out), dimension(:), allocatable :: dst
    include "include/copy_alloc_1d_impl.f90"
end subroutine

!-------------------------------------------------------------------------------
! Routines for 2d array

pure subroutine copy_alloc_2d_real32 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: PREC = real32
    integer, parameter :: ND = 2

    real (PREC), intent(in), dimension(:,:), optional :: src
    real (PREC), intent(out), dimension(:,:), allocatable :: dst
    include "include/copy_alloc_2d_impl.f90"
end subroutine

pure subroutine copy_alloc_2d_real64 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: PREC = real64
    integer, parameter :: ND = 2

    real (PREC), intent(in), dimension(:,:), optional :: src
    real (PREC), intent(out), dimension(:,:), allocatable :: dst
    include "include/copy_alloc_2d_impl.f90"
end subroutine

pure subroutine copy_alloc_2d_int32 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: INTSIZE = int32
    integer, parameter :: ND = 2

    integer (INTSIZE), intent(in), dimension(:,:), optional :: src
    integer (INTSIZE), intent(out), dimension(:,:), allocatable :: dst
    include "include/copy_alloc_2d_impl.f90"

end subroutine

pure subroutine copy_alloc_2d_int64 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: INTSIZE = int64
    integer, parameter :: ND = 2

    integer (INTSIZE), intent(in), dimension(:,:), optional :: src
    integer (INTSIZE), intent(out), dimension(:,:), allocatable :: dst
    include "include/copy_alloc_2d_impl.f90"

end subroutine

!-------------------------------------------------------------------------------
! Routines for 3d array

pure subroutine copy_alloc_3d_real32 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: PREC = real32
    integer, parameter :: ND = 3

    real (PREC), intent(in), dimension(:,:,:), optional :: src
    real (PREC), intent(out), dimension(:,:,:), allocatable :: dst
    include "include/copy_alloc_3d_impl.f90"
end subroutine

pure subroutine copy_alloc_3d_real64 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: PREC = real64
    integer, parameter :: ND = 3

    real (PREC), intent(in), dimension(:,:,:), optional :: src
    real (PREC), intent(out), dimension(:,:,:), allocatable :: dst
    include "include/copy_alloc_3d_impl.f90"
end subroutine

pure subroutine copy_alloc_3d_int32 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: INTSIZE = int32
    integer, parameter :: ND = 3

    integer (INTSIZE), intent(in), dimension(:,:,:), optional :: src
    integer (INTSIZE), intent(out), dimension(:,:,:), allocatable :: dst
    include "include/copy_alloc_3d_impl.f90"

end subroutine

pure subroutine copy_alloc_3d_int64 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: INTSIZE = int64
    integer, parameter :: ND = 3

    integer (INTSIZE), intent(in), dimension(:,:,:), optional :: src
    integer (INTSIZE), intent(out), dimension(:,:,:), allocatable :: dst
    include "include/copy_alloc_3d_impl.f90"

end subroutine
!-------------------------------------------------------------------------------
! Routines for 4d array

pure subroutine copy_alloc_4d_real32 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: PREC = real32
    integer, parameter :: ND = 4

    real (PREC), intent(in), dimension(:,:,:,:), optional :: src
    real (PREC), intent(out), dimension(:,:,:,:), allocatable :: dst
    include "include/copy_alloc_4d_impl.f90"
end subroutine

pure subroutine copy_alloc_4d_real64 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: PREC = real64
    integer, parameter :: ND = 4

    real (PREC), intent(in), dimension(:,:,:,:), optional :: src
    real (PREC), intent(out), dimension(:,:,:,:), allocatable :: dst
    include "include/copy_alloc_4d_impl.f90"
end subroutine

pure subroutine copy_alloc_4d_int32 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: INTSIZE = int32
    integer, parameter :: ND = 4

    integer (INTSIZE), intent(in), dimension(:,:,:,:), optional :: src
    integer (INTSIZE), intent(out), dimension(:,:,:,:), allocatable :: dst
    include "include/copy_alloc_4d_impl.f90"

end subroutine

pure subroutine copy_alloc_4d_int64 (src, dst)
    !*  COPY_ALLOC implements a routine similar to MOVE_ALLOC, but leaves
    !   the SRC argument unchanged. If the source and destination arrays have
    !   the same shape, the source data is copied directly into the destination
    !   array, instead of re-allocating the destination array.
    !
    !   Unlike MOVE_ALLOC, COPY_ALLOC does not (and cannot) modify any pointers
    !   to SRC.
    !
    !   If SRC is either missing (or not allocated, which is interpreted as
    !   not being present in Fortran 2008), DST becomes unallocated
    !   on exit.

    integer, parameter :: INTSIZE = int64
    integer, parameter :: ND = 4

    integer (INTSIZE), intent(in), dimension(:,:,:,:), optional :: src
    integer (INTSIZE), intent(out), dimension(:,:,:,:), allocatable :: dst
    include "include/copy_alloc_4d_impl.f90"

end subroutine



end module
