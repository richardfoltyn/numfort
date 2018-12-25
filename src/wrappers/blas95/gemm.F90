

subroutine __APPEND(gemm,__PREC) (a, b, c, transa, transb, alpha, beta)
    !*  GEMM provides a Fortran 95 wrapper for the Fortran 77 BLAS routines
    !   ?GEMM which is API-compatible with the wrapper included in
    !   Intel's MKL.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:), contiguous :: a
    real (PREC), intent(in), dimension(:,:), contiguous :: b
    real (PREC), intent(out), dimension(:,:), contiguous :: c
    character (1), intent(in), optional :: transa
    character (1), intent(in), optional :: transb
    real (PREC), intent(in), optional :: alpha
    real (PREC), intent(in), optional :: beta

    real (PREC) :: lalpha, lbeta
    character (1) :: ltransa, ltransb
    integer :: m, n, k, lda, ldb, ldc

    lalpha = 1.0_PREC
    lbeta = 0.0_PREC
    ltransa = 'N'
    ltransb = 'N'

    if (present(alpha)) lalpha = alpha
    if (present(beta)) lbeta = beta
    if (present(transa)) ltransa = transa
    if (present(transb)) ltransb = transb

    m = size(a, 1)
    if (ltransa == 'N' .or. ltransa == 'n') then
        m = size(a, 1)
        k = size(a, 2)
    else
        m = size(a, 2)
        k = size(a, 1)
    end if

    if (ltransb == 'N' .or. ltransb == 'n') then
        n = size(b, 2)
    else
        n = size(b, 1)
    end if

    lda = size(a, 1)
    ldb = size(b, 1)
    ldc = size(c, 1)

    call GEMM (ltransa, ltransb, m, n, k, lalpha, a, lda, b, ldb, lbeta, c, ldc)

end subroutine
