module lapack_interfaces

    implicit none
    private

    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: SP = kind(1.0e0)
    integer, parameter :: COMPLEX_SP = kind((1.0e0, 1.0e0))
    integer, parameter :: COMPLEX_DP = kind((1.0d0, 1.0d0))

    public :: GETRF
    public :: GETRI
    public :: GETRS

    interface GETRF
        procedure sgetrf, dgetrf, cgetrf, zgetrf
    end interface

    interface GETRI
        procedure SGETRI, DGETRI, CGETRI, ZGETRI
    end interface
    
    interface GETRS
        procedure SGETRS, DGETRS, CGETRS, ZGETRS
    end interface
    
    interface
        subroutine sgetrf (m, n, a, lda, ipiv, info)
            import SP
            integer :: m, n, lda, info
            real (SP) :: a(lda, *)
            integer :: ipiv(*)
        end subroutine

        subroutine dgetrf (m, n, a, lda, ipiv, info)
            import DP
            integer :: m, n, lda, info
            real (DP) :: a(lda, *)
            integer :: ipiv(*)
        end subroutine

        subroutine cgetrf (m, n, a, lda, ipiv, info)
            import COMPLEX_SP
            integer :: m, n, lda, info
            complex (COMPLEX_SP) :: a(lda, *)
            integer :: ipiv(*)
        end subroutine

        subroutine zgetrf (m, n, a, lda, ipiv, info)
            import COMPLEX_DP
            integer :: m, n, lda, info
            complex (COMPLEX_DP) :: a(lda, *)
            integer :: ipiv(*)
        end subroutine
    end interface

    interface
        subroutine SGETRI (n, a, lda, ipiv, work, lwork, info)
            import SP
            integer     :: info, lda, lwork, n
            integer     :: ipiv(*)
            real (SP)   :: a(lda, *), work(*)
        end subroutine

        subroutine DGETRI (n, a, lda, ipiv, work, lwork, info)
            import DP
            integer     :: info, lda, lwork, n
            integer     :: ipiv(*)
            real (DP)   :: a(lda, *), work(*)
        end subroutine

        subroutine CGETRI (n, a, lda, ipiv, work, lwork, info)
            import COMPLEX_SP
            integer             :: info, lda, lwork, n
            integer             :: ipiv(*)
            complex (COMPLEX_SP)   :: a(lda, *), work(*)
        end subroutine

        subroutine ZGETRI (n, a, lda, ipiv, work, lwork, info)
            import COMPLEX_DP
            integer             :: info, lda, lwork, n
            integer             :: ipiv(*)
            complex (COMPLEX_DP)   :: a(lda, *), work(*)
        end subroutine

    end interface


    interface 
        subroutine SGETRS (trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import SP
            character   :: trans
            integer     :: info, lda, ldb, n, nrhs
            integer     :: ipiv (*)
            real (SP)   :: a(lda,*), b(ldb,*)
        end subroutine
    
        subroutine DGETRS (trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import DP
            character   :: trans
            integer     :: info, lda, ldb, n, nrhs
            integer     :: ipiv (*)
            real (DP)   :: a(lda,*), b(ldb,*)
        end subroutine
    
        subroutine CGETRS (trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import COMPLEX_SP
            character               :: trans
            integer                 :: info, lda, ldb, n, nrhs
            integer                 :: ipiv (*)
            complex (COMPLEX_SP)    :: a(lda,*), b(ldb,*)
        end subroutine
    
        subroutine ZGETRS (trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import COMPLEX_DP
            character               :: trans
            integer                 :: info, lda, ldb, n, nrhs
            integer                 :: ipiv (*)
            complex (COMPLEX_DP)    :: a(lda,*), b(ldb,*)
        end subroutine
    
    end interface

end module
