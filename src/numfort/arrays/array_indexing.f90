module array_indexing

    use iso_fortran_env, only: int32, int64, real32, real64
    implicit none
    private 
    
    interface sub2ind
        module procedure sub2ind_int32, sub2ind_int64, sub2ind_1d_int32, sub2ind_1d_int64
    end interface
    
    interface ind2sub
        module procedure ind2sub_int32, ind2sub_int64, ind2sub_1d_int32, ind2sub_1d_int64
    end interface
    
    public :: sub2ind, ind2sub
    
contains
    
! ******************************************************************************
! SUB2IND

pure subroutine sub2ind_int64 (shp, sub_indices, lin_indices)
    integer, parameter :: INTSIZE = int64
    include "includes/sub2ind.f90"
end subroutine

pure subroutine sub2ind_int32 (shp, sub_indices, lin_indices)
    integer, parameter :: INTSIZE = int32
    include "includes/sub2ind.f90"
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
    include "includes/ind2sub.f90"    
end subroutine

pure subroutine ind2sub_int32 (shp, lin_indices, sub_indices)
    integer, parameter :: INTSIZE = int32
    include "includes/ind2sub.f90"
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

end module
    