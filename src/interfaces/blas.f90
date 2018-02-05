

module blas_interfaces
    
    implicit none
    private
    
    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: SP = kind(1.0e0)
    integer, parameter :: COMPLEX_SP = kind((1.0e0, 1.0e0))
    integer, parameter :: COMPLEX_DP = kind((1.0d0, 1.0d0))
    
    ! Exported BLAS 1 routines
    public :: AXPY
    public :: COPY
    public :: DOT
    public :: SCAL
    
    ! Exported BLAS 2 routines
    public :: GEMV
    public :: GER
    
    ! Exported BLAS 3 routines
    public :: GEMM
    
    interface AXPY
        procedure SAXPY, DAXPY, CAXPY, ZAXPY
    end interface

    interface COPY
        procedure SCOPY, DCOPY, CCOPY, ZCOPY
    end interface
    
    interface DOT
        procedure SDOT, DDOT
    end interface
    
    interface SCAL
        procedure SSCAL, DSCAL, CSCAL, ZSCAL
    end interface
    
    interface GEMV
        procedure SGEMV, DGEMV, CGEMV, ZGEMV
    end interface
    
    interface GER
        procedure SGER, DGER
    end interface
    
    interface GEMM
        procedure SGEMM, DGEMM, CGEMM, ZGEMM
    end interface

    interface
        subroutine SCOPY (n, x, incx, y, incy)
            import SP
            integer     :: n, incx, incy
            real (SP)   :: x(*), y(*)
        end subroutine

        subroutine DCOPY (n, x, incx, y, incy)
            import DP
            integer     :: n, incx, incy
            real (DP)   :: x(*), y(*)
        end subroutine

        subroutine CCOPY (n, x, incx, y, incy)
            import COMPLEX_SP
            integer                 :: n, incx, incy
            complex (COMPLEX_SP)    :: x(*), y(*)
        end subroutine

        subroutine ZCOPY (n, x, incx, y, incy)
            import COMPLEX_DP
            integer                 :: n, incx, incy
            complex (COMPLEX_DP)    :: x(*), y(*)
        end subroutine
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
    
   
    ! --------------------------------------------------------------------------
    ! BLAS 1: SCAL
    interface
        subroutine SSCAL (n, a, x, incx)
            import SP
            integer         :: n, incx
            real (SP)       :: a
            real (SP)       :: x(*)
        end subroutine
    
        subroutine DSCAL (n, a, x, incx)
            import DP
            integer         :: n, incx
            real (DP)       :: a
            real (DP)       :: x(*)
        end subroutine
    
        subroutine CSCAL (n, a, x, incx)
            import COMPLEX_SP
            integer                 :: n, incx
            complex (COMPLEX_SP)    :: a
            complex (COMPLEX_SP)    :: x(*)
        end subroutine
        
        subroutine CSSCAL (n, a, x, incx)
            import SP
            import COMPLEX_SP
            integer                 :: n, incx
            real (SP)               :: a
            complex (COMPLEX_SP)    :: x(*)
        end subroutine
    
        subroutine ZSCAL (n, a, x, incx)
            import COMPLEX_DP
            integer                 :: n, incx
            complex (COMPLEX_DP)    :: a
            complex (COMPLEX_DP)    :: x(*)
        end subroutine
        
        subroutine ZDSCAL (n, a, x, incx)
            import DP
            import COMPLEX_DP
            integer                 :: n, incx
            real (DP)               :: a
            complex (COMPLEX_DP)    :: x(*)
        end subroutine
    end interface
    
    ! -------------------------------------------------------------------------- 
    ! BLAS 1: AXPY
    
    interface
        subroutine SAXPY (n, a, x, incx, y, incy)
            import SP
            integer     :: n, incx, incy
            real (SP)   :: a
            real (SP)   :: x(*), y(*)
        end subroutine
    
        subroutine DAXPY (n, a, x, incx, y, incy)
            import DP
            integer     :: n, incx, incy
            real (DP)   :: a
            real (DP)   :: x(*), y(*)
        end subroutine
    
        subroutine CAXPY (n, a, x, incx, y, incy)
            import COMPLEX_SP
            integer                 :: n, incx, incy
            complex (COMPLEX_SP)    :: a
            complex (COMPLEX_SP)    :: x(*), y(*)
        end subroutine
    
        subroutine ZAXPY (n, a, x, incx, y, incy)
            import COMPLEX_DP
            integer                 :: n, incx, incy
            complex (COMPLEX_DP)    :: a
            complex (COMPLEX_DP)    :: x(*), y(*)
        end subroutine
    end interface
    
    ! -------------------------------------------------------------------------- 
    ! BLAS 2: GER
    
    interface 
    
        subroutine SGER (m, n, alpha, x, incx, y, incy, a, lda)
            import SP
            integer             :: m, n, incx, incy, lda
            real (SP)           :: alpha
            real (SP)           :: x(*), y(*), a(lda,*)
        end subroutine
    
        subroutine DGER (m, n, alpha, x, incx, y, incy, a, lda)
            import DP
            integer             :: m, n, incx, incy, lda
            real (DP)           :: alpha
            real (DP)           :: x(*), y(*), a(lda,*)
        end subroutine
    end interface
end module
