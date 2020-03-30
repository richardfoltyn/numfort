program minpack_hybr_examples
    !*  Example for root_hybrd and root_hybrj, wrappers around MINPACK's HYBRD
    !   and HYBRJ routines to solve nonlinear equations system with or without
    !   a user-provided Jacobian.

    use numfort_optimize, workspace => workspace_real64, &
        optim_result => optim_result_real64
    use iso_fortran_env

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()

    contains



subroutine example1 ()
    !*  Find root of function F: R^3->R^3

    integer, parameter :: n = 3
    real (PREC), dimension(n) :: x, fx
    type (optim_result) :: res
    type (workspace) :: work
    ! error flag for Jacobian diagnostics
    real (PREC), dimension(n) :: err

    real (PREC), parameter :: x0 = 5.0

    x = x0
    ! (one) root of fcn1_jac is located at x = [2,3,4]
    call root_hybrj (fcn1_jac, x, fx, work=work, res=res)
    call print_report (res)

    ! Compare analytical and numerical derivative using chkder
    x = x0
    call chkder (fcn1_jac, x, err)
    print '(t3, "Gradient diagnostics: ", *(f10.7, :, ","))', err

    ! Compare to results obtained from hybrd using numerical differentiation
    x = x0
    call root_hybrd (fcn1, x, fx, work=work, res=res)
    call print_report (res)

end subroutine



pure subroutine fcn1_jac (x, fx, jac, task)
    real (PREC), dimension(:), intent(in), contiguous :: x
    real (PREC), dimension(:), intent(out), contiguous :: fx
    real (PREC), dimension(:,:), intent(out), contiguous :: jac
    integer, intent(inout) :: task

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



pure subroutine fcn1 (x, fx)
    real (PREC), dimension(:), intent(in), contiguous :: x
    real (PREC), dimension(:), intent(out), contiguous :: fx

    real (PREC), dimension(size(x), size(x)) :: jac
    integer :: task

    task = 1

    call fcn1_jac (x, fx, jac, task)
end subroutine



subroutine print_report (res)
    type (optim_result), intent(in) :: res

    integer, save :: ii = 1

    print "('#', t3, 'Example ', i0)", ii
    print "(t3, 'Status code(s): ', a)", char(res%status)
    if (.not. res%success .and. len_trim(res%msg) > 0) then
        print "(t3, 'Message: ', a)", res%msg
    end if
    write (OUTPUT_UNIT, advance='no', &
        fmt="(t3, 'Solution:', t23, '[', t25, *(t25, 5(f12.8, :, ', '), :, /))") res%x
    print *, ']'
    write (OUTPUT_UNIT, advance='no', &
        fmt="(t3, 'Function value:', t23, '[', t25, *(t25, 5(f12.8, :, ', '), :, /))") res%fx
    print *, ']'
    print "(t3, 'Number of function evaluations: ', i0)", res%nfev

    ii = ii + 1

end subroutine

end
