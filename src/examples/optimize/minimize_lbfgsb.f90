program lbfgsb

    use numfort_optimize
    use iso_fortran_env

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()

contains


subroutine example1 ()
    type (workspace) :: ws
    type (optim_result) :: res
    integer, parameter :: n = 25, m = 5, iprint = OPTIM_PRINT_MINIMAL
    ! boundaries
    real (PREC), dimension(n) :: x, lbnd, ubnd

    ! starting point
    x = 3.0_PREC
    ! lower and upper bounds for odd indices
    lbnd(1:n:2) = 1.0_PREC
    ubnd(1:n:2) = 100_PREC
    lbnd(2:n:2) = -100_PREC
    ubnd(2:n:2) = 100_PREC

    call minimize_lbfgsb (fobj1, x, grad1, m=m, lbounds=lbnd, ubounds=ubnd, &
        iprint=iprint, work=ws, res=res)

    print "('Function value at minimum: ', en22.15e2)", res%fx_opt
    print "('Number of iterations: ', i0)", res%nit
    print "('Number of function evaluations: ', i0)", res%nfev
    write (OUTPUT_UNIT, advance='no', &
        fmt="('Optimum located at: [', t23, *(t23, 5(f6.4, :, ', '), :, /))") res%x_opt
    print *, ' ]'

end subroutine

! fobj1 corresponds to the objective function used in the example in
! driver1.f90 in the original L-BGFS-B package
pure function fobj1 (x) result(fx)
    real (PREC), intent(in), dimension(:) :: x
    real (PREC) :: fx

    integer :: i, n

    n = size(x)

    fx = .25d0*( x(1)-1.d0 )**2
    do i = 2, n
        fx = fx + (x(i)-x(i-1)**2)**2
    end do
    fx = 4 * fx
end function

! gradient of fobj1
pure subroutine grad1 (x, g)
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: g

    real (PREC) :: t1, t2
    integer :: i, n

    n = size(x)

    t1  = x(2) - x(1)**2
    g(1) = 2.d0*(x(1) - 1.d0) - 1.6d1*x(1)*t1
    do i = 2, n-1
       t2 = t1
       t1 = x(i+1) - x(i)**2
       g(i) = 8.d0*t2 - 1.6d1*x(i)*t1
    end do
    g(n) = 8.d0*t1

end subroutine



end program
