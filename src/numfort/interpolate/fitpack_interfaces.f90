! Module provides explicit interfaces for F77 fitpack routines.
! These should facility better compile-time error checking when calling
! these functions.

module fitpack_interfaces

    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter :: PREC = real64

    interface
        subroutine curfit_if (iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
            import PREC
            real (PREC) :: xb,xe,s,fp
            integer :: iopt,m,k,nest,n,lwrk,ier
            ! array arguments
            real (PREC) :: x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
            integer :: iwrk(nest)
        end subroutine

        subroutine splev_if (t,n,c,k,x,y,m,e,ier)
            import PREC
            integer :: n, k, m, e, ier
            real (PREC) :: t(n), c(n), x(m), y(m)
        end subroutine

        subroutine splder_if (t,n,c,k,nu,x,y,m,e,wrk,ier)
            import PREC
            integer :: n,k,nu,m,e,ier
            real (PREC) :: t(n),c(n),x(m),y(m),wrk(n)
        end subroutine

    end interface

    public :: curfit_if, splev_if, splder_if

contains

end module
