module numfort_arrays_create

    use, intrinsic :: iso_fortran_env
    use numfort_core
    implicit none
    private

    public :: arange, linspace, powerspace
    public :: diag, diag_matrix, identity
    public :: vander

    interface arange
        procedure arange_impl_int32, arange_impl_int64, arange_int32, arange_int64
    end interface

    interface linspace
        procedure linspace_real32, linspace_real64
    end interface

    interface powerspace
        procedure powerspace_real32, powerspace_real64, &
            powerspace_real32_int32, powerspace_real64_int32
    end interface

    interface diag
        procedure diag_real32, diag_real64
    end interface

    interface diag_matrix
        procedure diag_matrix_real32, diag_matrix_real64
    end interface

    interface identity
        procedure identity_real32, identity_real64
    end interface

    interface vander
        procedure vander_scalar_real32, vander_scalar_real64, &
            vander_real32, vander_real64
    end interface

contains

! ******************************************************************************
! LINSPACE procedures
pure subroutine linspace_real64(x, xfrom, xto, step, res_n, res_step)
    integer, parameter :: PREC = real64
    include "include/linspace_impl.f90"
end subroutine

pure subroutine linspace_real32(x, xfrom, xto, step, res_n, res_step)
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

! ------------------------------------------------------------------------------
! POWERSPACE routines

pure subroutine powerspace_real32 (x, xmin, xmax, pow)
    !*  POWERSPACE returns a sequence of points obtained by taking the power
    !   of a sequence of uniformly spaced points on [0,1] and applying
    !   an affine transformation to match the given start and end point.
    !
    !   More specifically, each element in X is computed as
    !       x(i) = xmin + (xmax-xmin) * u(i) ** pow
    !   where u(i) is the corresponding element on a uniformly-spaced
    !   sequence on [0,1].
    integer, parameter :: PREC = real32
    include "include/powerspace_impl.f90"
end subroutine

pure subroutine powerspace_real64 (x, xmin, xmax, pow)
    !*  POWERSPACE returns a sequence of points obtained by taking the power
    !   of a sequence of uniformly spaced points on [0,1] and applying
    !   an affine transformation to match the given start and end point.
    !
    !   More specifically, each element in X is computed as
    !       x(i) = xmin + (xmax-xmin) * u(i) ** pow
    !   where u(i) is the corresponding element on a uniformly-spaced
    !   sequence on [0,1].
    integer, parameter :: PREC = real64
    include "include/powerspace_impl.f90"
end subroutine

pure subroutine powerspace_real32_int32 (x, xmin, xmax, pow)
    !*  POWERSPACE returns a sequence of points obtained by taking the power
    !   of a sequence of uniformly spaced points on [0,1] and applying
    !   an affine transformation to match the given start and end point.
    !
    !   More specifically, each element in X is computed as
    !       x(i) = xmin + (xmax-xmin) * u(i) ** pow
    !   where u(i) is the corresponding element on a uniformly-spaced
    !   sequence on [0,1].
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:) :: x
        !*  Array to store "power-spaced" sequence of points
    real (PREC), intent(in) :: xmin
        !*  Starting value
    real (PREC), intent(in) :: xmax
        !*  Endpoint value
    integer, intent(in) :: pow
        !*  Exponent used to create power-spaced sequence

    call powerspace (x, xmin, xmax, real(pow, PREC))
end subroutine

pure subroutine powerspace_real64_int32 (x, xmin, xmax, pow)
    !*  POWERSPACE returns a sequence of points obtained by taking the power
    !   of a sequence of uniformly spaced points on [0,1] and applying
    !   an affine transformation to match the given start and end point.
    !
    !   More specifically, each element in X is computed as
    !       x(i) = xmin + (xmax-xmin) * u(i) ** pow
    !   where u(i) is the corresponding element on a uniformly-spaced
    !   sequence on [0,1].
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:) :: x
        !*  Array to store "power-spaced" sequence of points
    real (PREC), intent(in) :: xmin
        !*  Starting value
    real (PREC), intent(in) :: xmax
        !*  Endpoint value
    integer, intent(in) :: pow
        !*  Exponent used to create power-spaced sequence

    call powerspace (x, xmin, xmax, real(pow, PREC))
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
    include "include/vander_impl.f90"
end subroutine

pure subroutine vander_real64 (x, xp, increasing)
    !*  VANDER creates a Vandermonde matrix (or vector if the input is scalar)
    !   where columns of the output array are powers of input values, starting
    !   at 0. Power are by default arranged in decreasing order.

    integer, parameter :: PREC = real64
    include "include/vander_impl.f90"
end subroutine

end module
