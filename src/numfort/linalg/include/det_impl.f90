! Type/kind-independent implementation of matrix determinant function
! This uses LAPACK's LU factorization and computes the determinant as the
! product of diagonal elements of the resulting matrix, adjusting for
! the correct sign.
! Note: we need to create copy of input matrix as it otherwise would be
! overwritten by LAPACK's getrf()

real (PREC), intent (in), dimension(:,:) :: a
real (PREC) :: d

real (PREC), dimension(:,:), allocatable :: acopy
integer, dimension(:), allocatable :: ipiv

integer :: pow, n, linfo, i

if (size(a, 1) /= size(a, 2)) then
    write (ERROR_UNIT, *) "DET: a must be square matrix"
    return
end if

n = size(a, 1)

! create copy as a is overwritten by getrf
allocate (acopy(n,n), source=a)
allocate (ipiv(n))

ipiv = 0

call lapack_getrf(n, n, acopy, n, ipiv, linfo)

if (linfo /= 0) then
    write (ERROR_UNIT, *) "DET: Error while computing LU factorization"
    return
end if

pow = 0
d = 1.0_PREC

do i=1,n
    d = d * acopy(i,i)
    if (ipiv(i) /= i) pow = pow + 1
end do

! Correct sign
d = d * (-1) ** pow
