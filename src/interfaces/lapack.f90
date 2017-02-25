module lapack_interfaces

    implicit none
    private

    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: SP = kind(1.0e0)

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
            import SP
            integer :: m, n, lda, info
            complex (SP) :: a(lda, *)
            integer :: ipiv(*)
        end subroutine

        subroutine zgetrf (m, n, a, lda, ipiv, info)
            import DP
            integer :: m, n, lda, info
            complex (DP) :: a(lda, *)
            integer :: ipiv(*)
        end subroutine
    end interface

    interface lapack_getrf
        procedure sgetrf, dgetrf, cgetrf, zgetrf
    end interface

    public :: lapack_getrf

end module
