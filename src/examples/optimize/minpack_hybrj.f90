! Example for root_hybrj, a wrapper around MINPACK's hybrj routine to
! solve nonlinear equations system with a user-provided Jacobian.
program minpack_examples

    use numfort_common
    use numfort_optimize
    use iso_fortran_env

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()
contains

subroutine example1 ()

    integer, parameter :: n = 3
    real (PREC), dimension(n) :: x, fx
    type (optim_result) :: res
    type (workspace) :: work

    x = 1.0
    ! (one) root of func1 is located at x = [2,3,4]
    call root_hybrj (func1, x, fx, work=work, res=res)
    call print_report (res)

end subroutine

pure subroutine func1 (x, fx, jac, task)
    real (PREC), dimension(:), intent(in) :: x
    real (PREC), dimension(:), intent(out) :: fx
    real (PREC), dimension(:,:), intent(out) :: jac
    integer, intent(in out) :: task

    if (task == 1) then
        fx(1) = x(1)**2 + 2*x(2) - 2*x(3) - 2
        fx(2) = 3*x(1) - x(2)**2 + x(3) - 1
        fx(3) = 5*x(1) - 2*x(2) - x(3)**2 + 12
    else if (task == 2) then
        jac(1, :) = [2*x(1), 2.0d0, -2.0d0]
        jac(2, :) = [3.0d0, -2*x(2), 1.0d0]
        jac(3, :) = [5.0d0, -2.0d0, -2*x(3)]
    end if
end subroutine


subroutine print_report (res)
    class (optim_result), intent(in) :: res

    integer, save :: ii = 1

    print "('#', t3, 'Example ', i0)", ii
    print "(t3, 'Status code(s): ', a)", char(res%status)
    if (.not. res%success .and. len_trim(res%msg) > 0) then
        print "(t3, 'Message: ', a)", res%msg
    end if
    write (OUTPUT_UNIT, advance='no', &
        fmt="(t3, 'Solution:', t23, '[', t25, *(t25, 5(f6.4, :, ', '), :, /))") res%x
    print *, ']'
    write (OUTPUT_UNIT, advance='no', &
        fmt="(t3, 'Function value:', t23, '[', t25, *(t25, 5(f6.4, :, ', '), :, /))") res%fx
    print *, ']'
    print "(t3, 'Number of function evaluations: ', i0)", res%nfev

    ii = ii + 1

end subroutine

end
