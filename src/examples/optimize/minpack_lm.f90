program minpack_lm_demo
    !*  Example code using wrappers for MINPACKS LMDER and LMDIF
    !   routines for least-quares root finding.

    use numfort_optimize, workspace => workspace_real64, &
        optim_result => optim_result_real64
    use iso_fortran_env

    implicit none

    integer, parameter :: PREC = real64

    character (len=10) :: msg

    call example1 ()

contains

subroutine example1 ()

    ! Find root in Least-squares sense of a function F: R^2 -> R^4
    real (PREC), dimension(2) :: x
    real (PREC), dimension(4) :: fx

    type (optim_result) :: res
    type (workspace) :: work

    ! (one) solution for fcn1 is (2, 3)
    x = 0.0
    call root_lmdif (fcn1, x, fx, work=work, res=res)
    call print_report (res)

    x = 0.0
    call root_lmder (fcn1_jac, x, fx, work=work, res=res)
    call print_report (res)

end subroutine

pure subroutine fcn1_jac (x, fx, jac, task)
    real (PREC), dimension(:), intent(in) :: x
    real (PREC), dimension(:), intent(out) :: fx
    real (PREC), dimension(:,:), intent(out) :: jac
    integer, intent(in out) :: task

    real (PREC) :: x1, x2

    x1 = x(1)
    x2 = x(2)

    if (task == 1) then
        fx(1) = 2 * x1 - 2*x2 + 2
        fx(2) = x1**2 - 2*x1*x2 + 8
        fx(3) = x1 + 3*x2**2 - 29
        fx(4) = -2*x1 + x1*x2 - 2*x2 + 4
    else if (task == 2) then
        jac(1,:) = [2.0d0, -2.0d0]
        jac(2,:) = [2*x1-2*x2, -2*x1]
        jac(3,:) = [1.0d0, 6*x2]
        jac(4,:) = [-2+x2, x1-2]
    end if
end subroutine

pure subroutine fcn1 (x, fx)
    real (PREC), dimension(:), intent(in) :: x
    real (PREC), dimension(:), intent(out) :: fx

    real (PREC), dimension(0,0) :: jac
    integer :: task

    task = 1
    call fcn1_jac (x, fx, jac, task)

end subroutine


subroutine print_report (res)
    type (optim_result), intent(in) :: res

    integer, save :: ii = 1

    print "('#', t3, 'Example ', i0)", ii
    print "(t3, 'Status code(s): ', a)", char(res%status)
    if (len_trim(res%msg) > 0) then
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
