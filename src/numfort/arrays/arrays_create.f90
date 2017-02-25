module numfort_arrays_create

    use, intrinsic :: iso_fortran_env
    implicit none
    private

    interface arange
        module procedure arange_impl_int32, arange_impl_int64, arange_int32, arange_int64
    end interface

    interface linspace
        module procedure linspace_real32, linspace_real64
    end interface

    interface diag
        module procedure diag_real32, diag_real64
    end interface

    interface diag_matrix
        module procedure diag_matrix_real32, diag_matrix_real64
    end interface

    interface identity
        module procedure identity_real32, identity_real64
    end interface

    public :: arange, linspace
    public :: diag, diag_matrix, identity

contains

! ******************************************************************************
! LINSPACE procedures
subroutine linspace_real64(x, xfrom, xto)
    integer, parameter :: PREC = real64
    include "include/linspace_impl.f90"
end subroutine

subroutine linspace_real32(x, xfrom, xto)
    integer, parameter :: PREC = real32
    include "include/linspace_impl.f90"
end subroutine

! ******************************************************************************
! ARANGE procedures

pure subroutine arange_int32 (x)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:) :: x

    call arange (x, 1_INTSIZE, 1_INTSIZE)
end subroutine

pure subroutine arange_int64 (x)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:) :: x

    call arange (x, 1_INTSIZE, 1_INTSIZE)
end subroutine

pure subroutine arange_impl_int32 (x, ifrom, step)
    integer, parameter :: INTSIZE = int32
    include "include/arange_impl.f90"
end subroutine

pure subroutine arange_impl_int64 (x, ifrom, step)
    integer, parameter :: INTSIZE = int64
    include "include/arange_impl.f90"
end subroutine

! ******************************************************************************
! IDENTITY matrix creation

subroutine identity_real64(mat)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:) :: mat
    integer :: i

    mat = 0_PREC
    forall (i=1:size(mat, 1)) mat(i,i) = 1.0_PREC

end subroutine

subroutine identity_real32(mat)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:) :: mat
    integer :: i

    mat = 0_PREC
    forall (i=1:size(mat, 1)) mat(i,i) = 1.0_PREC

end subroutine

! ******************************************************************************
! DIAG - creating square matrices from diagonals
pure subroutine diag_matrix_real64(v, mat)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:) :: v
    real (PREC), intent(out), dimension(size(v), size(v)) :: mat

    integer :: i

    mat = 0.0_PREC
    forall (i=1:size(v)) mat(i,i) = v(i)
end subroutine

pure subroutine diag_matrix_real32(v, mat)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:) :: v
    real (PREC), intent(out), dimension(size(v), size(v)) :: mat

    integer :: i

    mat = 0.0_PREC
    forall (i=1:size(v)) mat(i,i) = v(i)

end subroutine

! ******************************************************************************
! DIAG - extracting diagonals from square matrices

pure subroutine diag_real64(mat, v)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:) :: mat
    real (PREC), intent(out), dimension(size(mat, 1)) :: v

    integer :: i

    forall (i=1:size(mat, 1)) v(i) = mat(i,i)

end subroutine

pure subroutine diag_real32(mat, v)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:) :: mat
    real (PREC), intent(out), dimension(size(mat, 1)) :: v

    integer :: i

    forall (i=1:size(mat, 1)) v(i) = mat(i,i)

end subroutine

end module
