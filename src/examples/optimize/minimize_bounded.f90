

program minimize_bounded_example

    use, intrinsic :: iso_fortran_env
    use numfort_optimize, optim_result => optim_result_real64

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()
    call example2 ()

    contains


subroutine example1 ()

    real (PREC) :: a, b, x, fx
    type (optim_result) :: res

    ! Test with interior minimium
    a = 0
    b = 4

    call minimize_bounded (fobj, a, b, x, res=res)
    fx = res%fx(1)
    print '(tr1, "Min. found at: ", f0.4, "; f(x)=", f0.4, "; true x_min: ", f0.4, "; nfev=", i0)', x, fx, 2.0, res%nfev

    ! Test with min. at boundary

    a = -5
    b = 1
    call minimize_bounded (fobj, a, b, x, res=res)
    fx = res%fx(1)
    print '(tr1, "Min. found at: ", f0.4, "; f(x)=", f0.4, "; true x_min: ", f0.4, "; nfev=", i0)', x, fx, 1.0, res%nfev
end subroutine


subroutine example2 ()

    real (PREC) :: a, b, x, fx
    type (optim_result) :: res
    real (PREC), dimension(2) :: args

    args(1) = 2.0
    args(2) = 1.0

    ! Test with interior minimium
    a = 0
    b = 4

    call minimize_bounded (fobj_args, a, b, x, args=args, res=res)
    fx = res%fx(1)
    print '(tr1, "Min. found at: ", f0.4, "; f(x)=", f0.4, "; true x_min: ", f0.4, "; nfev=", i0)', x, fx, 2.0, res%nfev

    ! Test with min. at boundary
    a = -5
    b = 1
    call minimize_bounded (fobj_args, a, b, x, args=args, res=res)
    fx = res%fx(1)
    print '(tr1, "Min. found at: ", f0.4, "; f(x)=", f0.4, "; true x_min: ", f0.4, "; nfev=", i0)', x, fx, 1.0, res%nfev
end subroutine


subroutine fobj (x, fx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx

    fx = (x - 2) ** 2
end subroutine

subroutine fobj_args (x, fx, args)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(in), dimension(:) :: args

    fx = (x - args(1)) ** 2 + args(2)
end subroutine

end program
