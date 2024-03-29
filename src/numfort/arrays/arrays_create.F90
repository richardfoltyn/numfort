module numfort_arrays_create

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks
    use numfort_core
    implicit none
    private

    public :: arange
    public :: diag, identity
    public :: vander
    public :: kron

    interface arange
        procedure arange_int32, arange_int64, arange_real32, arange_real64
    end interface

    interface diag
        procedure diag_real32, diag_real64
    end interface

    interface diag
        procedure diag_create_real32, diag_create_real64
    end interface

    interface identity
        procedure identity_int32, identity_int64, identity_real32, identity_real64
    end interface

    interface vander
        procedure vander_scalar_real32, vander_scalar_real64, &
            vander_real32, vander_real64
    end interface

    interface kron
        procedure kron_real32, kron_real64
    end interface

contains

! ******************************************************************************
! ARANGE procedures

pure subroutine arange_int32 (x, ifrom, step)
    integer, parameter :: INTSIZE = int32
#include "arange_impl.f90"
end subroutine

pure subroutine arange_int64 (x, ifrom, step)
    integer, parameter :: INTSIZE = int64
#include "arange_impl.f90"
end subroutine

pure subroutine arange_real32 (x, ifrom, step)
    integer, parameter :: PREC = real32
#include "arange_real_impl.f90"
end subroutine

pure subroutine arange_real64 (x, ifrom, step)
    integer, parameter :: PREC = real64
#include "arange_real_impl.f90"
end subroutine


! ******************************************************************************
! IDENTITY matrix creation

subroutine identity_int32 (mat)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:) :: mat

    integer :: i, n

    n = min(size(mat,1), size(mat,2))
    mat = 0_INTSIZE
    do i = 1, n
        mat(i,i) = 1_INTSIZE
    end do
end subroutine

subroutine identity_int64 (mat)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:,:) :: mat

    integer :: i, n

    n = min(size(mat,1), size(mat,2))
    mat = 0_INTSIZE
    do i = 1, n
        mat(i,i) = 1_INTSIZE
    end do
end subroutine

subroutine identity_real64(mat)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:) :: mat

    integer :: i, n

    n = min(size(mat,1), size(mat,2))
    mat = 0_PREC
    do i = 1, n
        mat(i,i) = 1.0_PREC
    end do
end subroutine

subroutine identity_real32(mat)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:) :: mat

    integer :: i, n

    n = min(size(mat,1), size(mat,2))
    mat = 0_PREC
    do i = 1, n
        mat(i,i) = 1.0_PREC
    end do
end subroutine

! ******************************************************************************
! DIAG - creating square matrices from diagonals


pure subroutine diag_create_real64 (v, res)
    !*  DIAG_CREATE creates a square matrix with a diagonal given by vector V.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:) :: v
    real (PREC), intent(out), dimension(:,:) :: res
        !*  Output array. If the size of MAT is not conformable with the
        !   diagonal matrix to be created, MAT is left unchanged.

    integer :: i, n
    integer :: shp(2)

    n = size(v)
    shp = n
    if (.not. has_shape (res, shp)) return

    res = 0.0_PREC
    do i = 1, n
        res(i,i) = v(i)
    end do
end subroutine

pure subroutine diag_create_real32 (v, res)
    !*  DIAG_CREATE creates a square matrix with a diagonal given by vector V.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:) :: v
    real (PREC), intent(out), dimension(:,:) :: res
        !*  Output array. If the size of MAT is not conformable with the
        !   diagonal matrix to be created, MAT is left unchanged.

    integer :: i, n
    integer :: shp(2)

    n = size(v)
    shp = n
    if (.not. has_shape (res, shp)) return

    res = 0.0_PREC
    do i = 1, n
        res(i,i) = v(i)
    end do

end subroutine

! ******************************************************************************
! DIAG - extracting diagonals from square matrices

pure subroutine diag_real64 (mat, res)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:) :: mat
    real (PREC), intent(out), dimension(:) :: res

    integer :: i, n

    n = min(size(mat,1), size(mat,2))
    if (size(res) < n) return

    do i = 1, n
        res(i) = mat(i,i)
    end do

end subroutine

pure subroutine diag_real32 (mat, res)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:) :: mat
    real (PREC), intent(out), dimension(:) :: res

    integer :: i, n

    n = min(size(mat,1), size(mat,2))
    if (size(res) < n) return

    do i = 1, n
        res(i) = mat(i,i)
    end do

end subroutine


!-------------------------------------------------------------------------------
! VANDER routines

pure subroutine vander_scalar_real32 (x, xp, increasing)
    !*  VANDER creates a Vandermonde matrix (or vector if the input is scalar)
    !   where columns of the output array are powers of input values, starting
    !   at 0. Power are by default arranged in decreasing order.

    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: x
        !   Scalar input value.
    real (PREC), intent(out), dimension(:), target :: xp
        !*  Vector of powers of X, where each column corresponds to an
        !   exponent in {0,1,...NP-1} where NP=SIZE(XP).
    logical, intent(in), optional :: increasing
        !*  If present and .TRUE., sort powers in increasing order, ie.
        !   first column contains X**0.0 = 1.0.

    real (PREC), dimension(:,:), pointer :: ptr_xp
    real (PREC), dimension(1) :: x1

    ptr_xp(1:1,1:size(xp)) => xp
    x1(1) = x
    call vander (x1, ptr_xp, increasing)
end subroutine

pure subroutine vander_scalar_real64 (x, xp, increasing)
    !*  VANDER creates a Vandermonde matrix (or vector if the input is scalar)
    !   where columns of the output array are powers of input values, starting
    !   at 0. Power are by default arranged in decreasing order.

    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: x
        !   Scalar input value.
    real (PREC), intent(out), dimension(:), target :: xp
        !*  Vector of powers of X, where each column corresponds to an
        !   exponent in {0,1,...NP-1} where NP=SIZE(XP).
    logical, intent(in), optional :: increasing
        !*  If present and .TRUE., sort powers in increasing order, ie.
        !   first column contains X**0.0 = 1.0.

    real (PREC), dimension(:,:), pointer :: ptr_xp
    real (PREC), dimension(1) :: x1

    ptr_xp(1:1,1:size(xp)) => xp
    x1(1) = x
    call vander (x1, ptr_xp, increasing)
end subroutine

pure subroutine vander_real32 (x, xp, increasing)
    !*  VANDER creates a Vandermonde matrix (or vector if the input is scalar)
    !   where columns of the output array are powers of input values, starting
    !   at 0. Power are by default arranged in decreasing order.

    integer, parameter :: PREC = real32
#include "vander_impl.f90"
end subroutine

pure subroutine vander_real64 (x, xp, increasing)
    !*  VANDER creates a Vandermonde matrix (or vector if the input is scalar)
    !   where columns of the output array are powers of input values, starting
    !   at 0. Power are by default arranged in decreasing order.

    integer, parameter :: PREC = real64
#include "vander_impl.f90"
end subroutine

!-------------------------------------------------------------------------------
! Kronecker product

pure subroutine kron_real32 (mat1, mat2, res)
    !*  KRON returns the Kronecker product of matrices MAT1 and MAT2.
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:) :: mat1
    real (PREC), intent(in), dimension(:,:) :: mat2
    real (PREC), intent(out), dimension(:,:), contiguous, target :: res
        !*  Output array. If the array size is non-conformable with
        !   the required dimensions to store the Kronecker product of
        !   MAT1 and MAT2, the routine leaves RES unchanged.

#include "kron_impl.F90"
end subroutine

pure subroutine kron_real64 (mat1, mat2, res)
    !*  KRON returns the Kronecker product of matrices MAT1 and MAT2.
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:) :: mat1
    real (PREC), intent(in), dimension(:,:) :: mat2
    real (PREC), intent(out), dimension(:,:), contiguous, target :: res
        !*  Output array. If the array size is non-conformable with
        !   the required dimensions to store the Kronecker product of
        !   MAT1 and MAT2, the routine leaves RES unchanged.

#include "kron_impl.F90"
end subroutine


end module
