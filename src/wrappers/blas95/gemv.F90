

subroutine __APPEND(gemv,__PREC) (a, x, y, alpha, beta, trans)
    !*  GEMV provides a Fortran 95 wrapper for the Fortran 77 BLAS routines
    !   ?GEMV which is API-compatible with the wrapper included in
    !   Intel's MKL.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:), contiguous :: a
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(inout), dimension(:), contiguous :: y
    real (PREC), intent(in), optional :: alpha
    real (PREC), intent(in), optional :: beta
    character (1), intent(in), optional :: trans

    character (1) :: ltrans
    real (PREC) :: lalpha, lbeta
    integer :: m, n, lda
    integer, parameter :: incx = 1, incy = 1

    lalpha = 1.0_PREC
    lbeta = 0.0_PREC
    ltrans = 'N'

    if (present(alpha)) lalpha = alpha
    if (present(beta)) lbeta = beta
    if (present(trans)) ltrans = trans

    m = size(a, 1)
    n = size(a, 2)
    lda = m

    call GEMV (ltrans, m, n, lalpha, a, lda, x, incx, lbeta, y, incy)

end subroutine
