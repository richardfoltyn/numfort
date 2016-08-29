! implementation of Matrix inversion subroutine

real (PREC), dimension(:,:), contiguous :: A
real (PREC), dimension(:,:), contiguous :: Ainv
integer :: info

intent(in) :: A
intent(out) :: Ainv, info
optional :: info

real(PREC), dimension(:), allocatable :: work  ! work array for LAPACK
integer, dimension(:), allocatable :: ipiv   ! pivot indices
integer :: n, linfo

! External procedures defined in LAPACK
external DGETRF
external DGETRI

! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
n = size(A,1)

allocate (work(n))
allocate (ipiv(n))

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
call DGETRF(n, n, Ainv, n, ipiv, linfo)

if (linfo /= 0) then
   write (ERROR_UNIT, *) 'Matrix is numerically singular!'
   goto 100
end if

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
call DGETRI(n, Ainv, n, ipiv, work, n, linfo)

if (linfo /= 0) then
   write (ERROR_UNIT, *) 'Matrix inversion failed!'
   goto 100
end if

100 if (present(info)) info = linfo
deallocate (work, ipiv)
