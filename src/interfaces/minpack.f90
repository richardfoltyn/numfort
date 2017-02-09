module minpack_interfaces

    use iso_fortran_env, only: real64
    implicit none

    integer, private, parameter :: PREC = real64

    interface
        subroutine hybrd (fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
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

        subroutine hybrj (fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode, &
                        factor,nprint,info,nfev,njev,r,lr,qtf,wa1,wa2, &
                        wa3,wa4)
            import PREC
            integer, intent(in) :: n, ldfjac, maxfev, mode, nprint, lr
            integer, intent(out) :: info, nfev, njev
            real (PREC), intent(in) :: xtol, factor
            real (PREC), intent(in out) :: x(n), fvec(n), fjac(ldfjac,n), &
                diag(n), r(lr), qtf(n),  wa1(n), wa2(n), wa3(n), wa4(n)
            external fcn
        end subroutine

        subroutine lmdif (fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn, &
                        diag,mode,factor,nprint,info,nfev,fjac,ldfjac, &
                        ipvt,qtf,wa1,wa2,wa3,wa4)

            import PREC
            external fcn
            integer, intent(in) :: m, n, maxfev, mode, nprint, ldfjac
            integer, intent(out) :: info, nfev
            integer, intent(out) :: ipvt(n)
            real (PREC), intent(in) :: ftol, xtol, gtol, epsfcn, factor
            real (PREC), intent(in) :: diag(n)
            real (PREC), intent(in out) :: x(n), fvec(n), fjac(ldfjac, n), &
                qtf(n), wa1(n), wa2(n), wa3(n), wa4(n)
        end subroutine

        subroutine lmder (fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, &
                        maxfev,diag,mode,factor,nprint,info,nfev,njev, &
                        ipvt,qtf,wa1,wa2,wa3,wa4)
            import PREC
            external fcn
            integer, intent(in) :: m, n, maxfev, mode, nprint, ldfjac
            integer, intent(out) :: info, nfev, njev
            integer, intent(out) :: ipvt(n)
            real (PREC), intent(in) :: ftol, xtol, gtol, factor
            real (PREC), intent(in) :: diag(n)
            real (PREC), intent(in out) :: x(n), fvec(n), fjac(ldfjac, n), &
                qtf(n), wa1(n), wa2(n), wa3(n), wa4(n)
        end subroutine
    end interface

end module
