module lapack_interfaces

    implicit none
    private

    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: SP = kind(1.0e0)
    integer, parameter :: COMPLEX_SP = kind((1.0e0, 1.0e0))
    integer, parameter :: COMPLEX_DP = kind((1.0d0, 1.0d0))

    public :: GELSD
    public :: GETRF
    public :: GETRI
    public :: GETRS
    public :: GESDD
    public :: GESV
    public :: GESVD
    public :: SYEVD

    interface GETRF
        procedure sgetrf, dgetrf, cgetrf, zgetrf
    end interface

    interface GETRI
        procedure SGETRI, DGETRI, CGETRI, ZGETRI
    end interface
    
    interface GETRS
        procedure SGETRS, DGETRS, CGETRS, ZGETRS
    end interface

    interface GESV
        !*  DGESV computes the solution to system of linear equations
        !   A * X = B for GE matrices.
        procedure SGESV, DGESV, CGESV, ZGESV, DSGESV, ZCGESV
    end interface

    interface GESVD
        !*  DGESVD computes the singular value decomposition (SVD) for
        !   GE matrices.
        procedure SGESVD, DGESVD, CGESVD, ZGESVD
    end interface

    interface GESDD
        !*  Computes the singular value decomposition of a general rectangular
        !   matrix using a divide and conquer method.
        procedure SGESDD, DGESDD, CGESDD, ZGESDD
    end interface GESDD

    interface GELSD
        procedure SGELSD, DGELSD, CGELSD, ZGELSD
    end interface

    interface SYEVD
        procedure SSYEVD, DSYEVD
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


    interface
        subroutine SGESV (n, nrhs, a, lda, ipiv, b, ldb, info)
            import SP
            integer     :: info, lda, ldb, n, nrhs
            integer     :: ipiv(*)
            real (SP)   :: a(lda,*), b(ldb,*)
        end subroutine

        subroutine DGESV (n, nrhs, a, lda, ipiv, b, ldb, info)
            import DP
            integer     :: info, lda, ldb, n, nrhs
            integer     :: ipiv(*)
            real (DP)   :: a(lda,*), b(ldb,*)
        end subroutine

        subroutine CGESV (n, nrhs, a, lda, ipiv, b, ldb, info)
            import COMPLEX_SP
            integer                 :: info, lda, ldb, n, nrhs
            integer                 :: ipiv(*)
            complex (COMPLEX_SP)    :: a(lda,*), b(ldb,*)
        end subroutine

        subroutine ZGESV (n, nrhs, a, lda, ipiv, b, ldb, info)
            import COMPLEX_DP
            integer                 :: info, lda, ldb, n, nrhs
            integer                 :: ipiv(*)
            complex (COMPLEX_DP)    :: a(lda,*), b(ldb,*)
        end subroutine

        subroutine DSGESV (n, nrhs, a, lda, ipiv, b, ldb, x, ldx, work, swork, iter, info)
            import DP
            import SP
            integer     :: info, iter, lda, ldb, ldx, n, nrhs
            integer     :: ipiv(*)
            real (SP)   :: swork(*)
            real (DP)   :: a(lda,*), b(ldb,*), work(n,*), x(ldx,*)
        end subroutine

        subroutine ZCGESV (n, nrhs, a, lda, ipiv, b, ldb, x, ldx, work, swork, rwork, iter, info)
            import DP
            import COMPLEX_SP
            import COMPLEX_DP
            integer                 :: info, iter, lda, ldb, ldx, n, nrhs
            integer                 :: ipiv(*)
            real (DP)               :: rwork(*)
            complex (COMPLEX_SP)    :: swork(*)
            complex (COMPLEX_DP)    :: a(lda,*), b(ldb,*), work(n,*), x(ldx,*)
        end subroutine
    end interface


    interface
        subroutine SGESVD (jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            import
            character               :: jobu, jobvt
            integer                 :: info, lda, ldu, ldvt, lwork, m, n
            real (SP)               :: a(lda, *), s(*), u(ldu, *), vt(ldvt, *), work(*)
        end subroutine

        subroutine DGESVD (jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            import
            character               :: jobu, jobvt
            integer                 :: info, lda, ldu, ldvt, lwork, m, n
            real (DP)               :: a(lda, *), s(*), u(ldu, *), vt(ldvt, *), work(*)
        end subroutine

        subroutine CGESVD (jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
            import
            character               :: jobu, jobvt
            integer                 :: info, lda, ldu, ldvt, lwork, m, n
            real (SP)               :: rwork(*), s(*)
            complex (COMPLEX_SP)    :: a(lda, *), u(ldu, *), vt(ldvt, *), work(*)
        end subroutine

        subroutine ZGESVD (jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
            import
            character               :: jobu, jobvt
            integer                 :: info, lda, ldu, ldvt, lwork, m, n
            real (DP)               :: rwork(*), s(*)
            complex (COMPLEX_DP)    :: a(lda, *), u(ldu, *), vt(ldvt, *), work(*)
        end subroutine
    end interface

    interface
        subroutine SGESDD (jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
            import
            character               :: jobz
            integer                 :: info, lda, ldu, ldvt, lwork, m, n
            integer                 :: iwork(*)
            real (SP)               :: a(lda, *), s(*), u(ldu, *), vt(ldvt, *), work(*)
        end subroutine

        subroutine DGESDD (jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
            import
            character               :: jobz
            integer                 :: info, lda, ldu, ldvt, lwork, m, n
            integer                 :: iwork(*)
            real (DP)               :: a(lda, *), s(*), u(ldu, *), vt(ldvt, *), work(*)
        end subroutine

        subroutine CGESDD (jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)
            import
            character               :: jobz
            integer                 :: info, lda, ldu, ldvt, lwork, m, n
            integer                 :: iwork(*)
            real (SP)               :: rwork(*), s(*)
            complex (COMPLEX_SP)    :: a(lda, *), u(ldu, *), vt(ldvt, *), worK(*)
        end subroutine

        subroutine ZGESDD (jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)
            import
            character               :: jobz
            integer                 :: info, lda, ldu, ldvt, lwork, m, n
            integer                 :: iwork(*)
            real (DP)               :: rwork(*), s(*)
            complex (COMPLEX_DP)    :: a(lda, *), u(ldu, *), vt(ldvt, *), work(*)
        end subroutine

    end interface

    interface
        subroutine SGELSD (m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
            import SP
            integer     :: info, lda, ldb, lwork, m, n, nrhs, rank
            real (SP)   :: rcond
            integer     :: iwork(*)
            real (SP)   :: a(lda,*), b(ldb,*), s(*), work(*)
        end subroutine

        subroutine DGELSD (m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
            import DP
            integer     :: info, lda, ldb, lwork, m, n, nrhs, rank
            real (DP)   :: rcond
            integer     :: iwork(*)
            real (DP)   :: a(lda,*), b(ldb,*), s(*), work(*)
        end subroutine

        subroutine CGELSD (m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info)
            import COMPLEX_SP
            import SP
            integer                 :: info, lda, ldb, lwork, m, n, nrhs, rank
            real (SP)               :: rcond
            integer                 :: iwork(*)
            real (SP)               :: rwork(*), s(*)
            complex (COMPLEX_SP)    :: a(lda,*), b(ldb,*), work(*)
        end subroutine

        subroutine ZGELSD (m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info)
            import COMPLEX_DP
            import DP
            integer                 :: info, lda, ldb, lwork, m, n, nrhs, rank
            real (DP)               :: rcond
            integer                 :: iwork(*)
            real (DP)               :: rwork(*), s(*)
            complex (COMPLEX_DP)    :: a(lda,*), b(ldb,*), work(*)
        end subroutine
    end interface


    interface
        subroutine SSYEVD (jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
            import
            integer                 :: lda, liwork, lwork, n, info
            character               :: jobz, uplo
            integer                 :: iwork(*)
            real (SP)               :: a(lda,*), w(*), work(*)
        end subroutine

        subroutine DSYEVD (jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
            import
            integer                 :: lda, liwork, lwork, n, info
            character               :: jobz, uplo
            integer                 :: iwork(*)
            real (DP)               :: a(lda,*), w(*), work(*)
        end subroutine
    end interface

end module
