program simplex_demo

    use numfort, only: pi
    use numfort_optimize, optim_result => optim_result_real64
    use iso_fortran_env

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()

    contains

! replicates example from t_minim.f90 that originally came with the simplex
! code.
subroutine example1 ()

    type (optim_result) :: res
    integer, parameter :: n = 3
    real (PREC) :: x(n), tol

    integer :: maxfun
    integer (NF_ENUM_KIND) :: iprint
    logical :: quad

    x = [real (PREC) :: 0, 1, 2]
    iprint = NF_PRINT_ALL
    maxfun = 250
    quad = .true.
    tol = 1d-04

    call minimize_simplex (fobj1, x, tol, maxfun, quad=quad, iprint=iprint, res=res)

    print "('Status codes(s): ', a, ' - ', a100)", char(res%status), res%msg
    print "('Function value at minimum: ', en22.15e2)", res%fx(1)
    write (OUTPUT_UNIT, advance='no', &
        fmt="('Optimum located at: [', t23, *(t23, 5(f6.4, :, ', '), :, /))") res%x
    print *, ' ]'

end subroutine

subroutine fobj1 (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx

    fx = - 1.0_PREC/(1.0_PREC + (x(1)-x(2)) ** 2) &
        - sin(pi * x(2) * x(3) * 0.5_PREC) &
        - exp(-((x(1)+x(3))/x(2) - 2.0_PREC) ** 2)
end subroutine

end program
