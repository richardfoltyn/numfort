module numfort_arrays_index
    !*  Collection of routines performing operations on array indices.

    use, intrinsic :: iso_fortran_env
    use numfort_arrays_create
    
    implicit none
    private

    public :: ind2sub
    public :: sub2ind
    public :: shape2sub

    !>  Converts an array of coordinate tuples into an array of flat indices.
    !   Inverse operation of `ind2sub()`.
    interface sub2ind
        module procedure sub2ind_int32, sub2ind_int64, sub2ind_1d_int32, sub2ind_1d_int64
    end interface

    !>  Convert an array of flat indices into an array of coordinate tuples.
    !   Inverse operation of `sub2ind()`.
    interface ind2sub
        module procedure ind2sub_int32, ind2sub_int64, ind2sub_1d_int32, ind2sub_1d_int64
    end interface

    !>  Create an array of coordinate tuples containing all indices for an array
    !   of given shape.
    interface shape2sub
        module procedure shape2sub_int32, shape2sub_int64
    end interface


contains

! ******************************************************************************
! SUB2IND

pure subroutine sub2ind_int64 (shp, sub_indices, lin_indices)
    integer, parameter :: INTSIZE = int64
#include "include/sub2ind_impl.F90"
end subroutine

pure subroutine sub2ind_int32 (shp, sub_indices, lin_indices)
    integer, parameter :: INTSIZE = int32
#include "include/sub2ind_impl.F90"
end subroutine

pure subroutine sub2ind_1d_int32 (shp, sub, ind)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:) :: shp
    integer (INTSIZE), intent(in), dimension(:) :: sub
    integer (INTSIZE), intent(out), dimension(:) :: ind

    ind = sub
end subroutine

pure subroutine sub2ind_1d_int64 (shp, sub, ind)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:) :: shp
    integer (INTSIZE), intent(in), dimension(:) :: sub
    integer (INTSIZE), intent(out), dimension(:) :: ind

    ind = sub
end subroutine

! ******************************************************************************
! IND2SUB

pure subroutine ind2sub_int64 (shp, lin_indices, sub_indices)
    integer, parameter :: INTSIZE = int64
#include "include/ind2sub_impl.F90"
end subroutine

pure subroutine ind2sub_int32 (shp, lin_indices, sub_indices)
    integer, parameter :: INTSIZE = int32
#include "include/ind2sub_impl.F90"
end subroutine

pure subroutine ind2sub_1d_int32 (shp, ind, sub)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:) :: shp
    integer (INTSIZE), intent(in), dimension(:) :: ind
    integer (INTSIZE), intent(out), dimension(:) :: sub

    sub = ind
end subroutine

pure subroutine ind2sub_1d_int64 (shp, ind, sub)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:) :: shp
    integer (INTSIZE), intent(in), dimension(:) :: ind
    integer (INTSIZE), intent(out), dimension(:) :: sub

    sub = ind
end subroutine

! ******************************************************************************
! SHAPE TO SUBINDEX

pure subroutine shape2sub_int64 (shp, sub)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:) :: shp
    integer (INTSIZE), intent(out), dimension(:, :) :: sub

    integer (INTSIZE), dimension(:), allocatable :: lidx

    allocate (lidx(product(shp)))

    call arange (lidx)
    call ind2sub (shp, lidx, sub)

    deallocate (lidx)

end subroutine

pure subroutine shape2sub_int32 (shp, sub)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:) :: shp
    integer (INTSIZE), intent(out), dimension(:, :) :: sub

    integer (INTSIZE), dimension(:), allocatable :: lidx

    allocate (lidx(product(shp)))

    call arange (lidx)
    call ind2sub (shp, lidx, sub)

    deallocate (lidx)

end subroutine
end module
