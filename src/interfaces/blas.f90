

module blas_interfaces
    
    implicit none
    private
    
    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: SP = kind(1.0e0)
    integer, parameter :: COMPLEX_SP = kind((1.0e0, 1.0e0))
    integer, parameter :: COMPLEX_DP = kind((1.0d0, 1.0d0))
    
    ! Exported BLAS 1 routines
    public :: DOT
    
    ! Exported BLAS 2 routines
    public :: GEMV
    
    ! Exported BLAS 3 routines
    public :: GEMM
    
    interface DOT
        procedure SDOT, DDOT
    end interface
    
    interface GEMM
        procedure SGEMM, DGEMM, CGEMM, ZGEMM
    end interface
    
    interface GEMV
        procedure SGEMV, DGEMV, CGEMV, ZGEMV
    end interface
    
    interface
        subroutine SGEMM (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import SP
            REAL (SP) ::    alpha,beta
            INTEGER   ::    k,lda,ldb,ldc,m,n
            CHARACTER ::    transa,transb
            REAL (SP) ::    a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
        
        subroutine DGEMM (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import DP
            REAL (DP) ::    alpha,beta
            INTEGER   ::    k,lda,ldb,ldc,m,n
            CHARACTER ::    transa,transb
            REAL (DP) ::    a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
        
        subroutine CGEMM (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import COMPLEX_SP
            COMPLEX (COMPLEX_SP) ::     alpha,beta
            INTEGER   ::                k,lda,ldb,ldc,m,n
            CHARACTER ::                transa,transb
            COMPLEX (COMPLEX_SP) ::     a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
        
        subroutine ZGEMM (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import COMPLEX_DP
            COMPLEX (COMPLEX_DP) ::     alpha,beta
            INTEGER   ::                k,lda,ldb,ldc,m,n
            CHARACTER ::                transa,transb
            COMPLEX (COMPLEX_DP) ::     a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
    end interface
    
    interface 
        subroutine SGEMV (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import SP
            character       :: trans
            real (SP)       :: alpha, beta
            integer         :: incx, incy, lda, m, n
            real (SP)       :: a(lda,*), x(*), y(*)
        end subroutine
    
        subroutine DGEMV (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import DP
            character       :: trans
            real (DP)       :: alpha, beta
            integer         :: incx, incy, lda, m, n
            real (DP)       :: a(lda,*), x(*), y(*)
        end subroutine
    
        subroutine CGEMV (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import SP
            character       :: trans
            complex (SP)    :: alpha, beta
            integer         :: incx, incy, lda, m, n
            complex (SP)    :: a(lda,*), x(*), y(*)
        end subroutine
    
        subroutine ZGEMV (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import DP
            character       :: trans
            complex (DP)    :: alpha, beta
            integer         :: incx, incy, lda, m, n
            complex (DP)    :: a(lda,*), x(*), y(*)
        end subroutine
    end interface
    
    interface
        function SDOT (n, x, incx, y, incy) result(res)
            import SP
            integer         :: n, incx, incy
            real (SP)       :: x(*), y(*)
            real (SP)       :: res
        end function
        
        function DDOT (n, x, incx, y, incy) result(res)
            import DP
            integer         :: n, incx, incy
            real (DP)       :: x(*), y(*)
            real (DP)       :: res
        end function
    end interface
end module
