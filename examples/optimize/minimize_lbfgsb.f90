program lbfgsb

    use numfort_optimize, workspace => workspace_real64, &
        optim_result => optim_result_real64
    use iso_fortran_env

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()
    call example2 ()

contains


subroutine example1 ()
    type (workspace) :: ws
    type (optim_result) :: res
    integer, parameter :: n = 25, m = 5
    integer (NF_ENUM_KIND), parameter :: iprint = NF_PRINT_MINIMAL
    ! boundaries
    real (PREC), dimension(n) :: x, lbnd, ubnd

    ! starting point
    x = 3.0_PREC
    ! lower and upper bounds for odd indices
    lbnd(1:n:2) = 1.0_PREC
    ubnd(1:n:2) = 100_PREC
    lbnd(2:n:2) = -100_PREC
    ubnd(2:n:2) = 100_PREC

    call minimize_lbfgsb (fobj1, grad1, x, m=m, lbounds=lbnd, ubounds=ubnd, &
        iprint=iprint, work=ws, res=res)
    call print_report (res)
end subroutine

! fobj1 corresponds to the objective function used in the example in
! driver1.f90 in the original L-BGFS-B package
pure subroutine fobj1 (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx

    integer :: i, n

    n = size(x)

    fx = .25d0*( x(1)-1.d0 )**2
    do i = 2, n
        fx = fx + (x(i)-x(i-1)**2)**2
    end do
    fx = 4 * fx
end subroutine

! gradient of fobj1
pure subroutine grad1 (x, g)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: g

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

! Identical to example 1, except that now be use a function argument
! that returns both the function value and the gradient.
subroutine example2 ()
    type (workspace) :: ws
    type (optim_result) :: res
    integer, parameter :: n = 25, m = 5
    integer (NF_ENUM_KIND), parameter :: iprint = NF_PRINT_MINIMAL
    ! boundaries
    real (PREC), dimension(n) :: x, lbnd, ubnd

    ! starting point
    x = 3.0_PREC
    ! lower and upper bounds for odd indices
    lbnd(1:n:2) = 1.0_PREC
    ubnd(1:n:2) = 100_PREC
    lbnd(2:n:2) = -100_PREC
    ubnd(2:n:2) = 100_PREC

    call minimize_lbfgsb (fobj_grad, x, m=m, lbounds=lbnd, ubounds=ubnd, &
        iprint=iprint, work=ws, res=res)

    call print_report (res)
end subroutine

subroutine fobj_grad (x, fx, g)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(out), dimension(:), contiguous ::  g

    call fobj1 (x, fx)
    call grad1 (x, g)
end subroutine

subroutine print_report(res)

    type (optim_result), intent(in) :: res
    integer, save :: ii = 1

    print "('#', t3, 'Example ', i0)", ii
    print "(t3, 'Exit status: ', a)", char(res%status)
    print "(t3, 'Function value at minimum: ', es10.2e3)", res%fx(1)
    print "(t3, 'Number of iterations: ', i0)", res%nit
    print "(t3, 'Number of function evaluations: ', i0)", res%nfev
    write (OUTPUT_UNIT, advance='no', &
        fmt="(t3, 'Optimum located at: [', t25, *(t25, 5(f6.4, :, ', '), :, /))") res%x
    print *, ' ]'

    ii = ii + 1

end subroutine



end program
