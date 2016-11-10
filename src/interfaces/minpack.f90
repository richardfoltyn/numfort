module minpack_interfaces

    use iso_fortran_env, only: real64
    implicit none

    integer, private, parameter :: PREC = real64

    interface
        subroutine hybrd_if (fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
                        mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr, &
                        qtf,wa1,wa2,wa3,wa4)
            import PREC
            integer, intent(in) :: n, maxfev, ml, mu, mode, nprint, ldfjac, lr
            integer, intent(out) :: info, nfev
            real (PREC), intent(in) :: xtol, epsfcn, factor
            real (PREC), intent(in out) :: x(n), fvec(n), diag(n), &
                fjac(ldfjac,n), r(lr), qtf(n), wa1(n), wa2(n), wa3(n), wa4(n)
            external fcn
        end subroutine

        subroutine hybrd1_if (fcn,n,x,fvec,tol,info,wa,lwa)
            import PREC
            integer, intent(in) :: n, lwa
            integer, intent(out) :: info
            real (PREC), intent(in) :: tol
            real (PREC), intent(in out) :: x(n), fvec(n), wa(lwa)
            external fcn
        end subroutine
    end interface

end module
