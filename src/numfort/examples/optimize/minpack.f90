program minpack_examples

    use numfort_optimize
    use iso_fortran_env

    implicit none

    integer, parameter :: PREC = real64

    character (len=10) :: msg
    read (*, *) msg

    call example1 ()
    call example2 ()

contains


! Examples for hybr() implemenetation of Powell's hybrid method
subroutine example1 ()

    real (PREC), dimension(2) :: x, fx
    type (optim_result) :: res
    type (workspace) :: work

    ! solution for func1 is x1 = 22/23 and x2 = 18/23
    x = 0.0

    call root_hybr (func1, x, fx, work=work, res=res)
    call print_report (res)

    ! nonlinear equation system
    ! (one?) solution for func2 == 0 is x1 = 1.21741 and x2 = 2.69888
    x = 2.0
    ! test call without passing in workspace object
    call root_hybr (func2, x, fx, res=res)
    call print_report (res)

end subroutine

! Examples for lstsq using the Levenberg-Marquardt algorithm
subroutine example2 ()

    real (PREC), dimension(2) :: x, fx
    type (optim_result) :: res
    type (workspace) :: work

    ! solution for func1 is x1 = 22/23 and x2 = 18/23
    x = 0.0

    call root_lstsq (func1, x, fx, work=work, res=res)
    call print_report (res)

    ! nonlinear equation system
    ! (one?) solution for func2 == 0 is x1 = 1.21741 and x2 = 2.69888
    x = 2.0
    ! test call without passing in workspace object
    call root_lstsq (func2, x, fx, res=res)
    call print_report (res)

end subroutine

pure subroutine func1 (x, fx)
    real (PREC), dimension(:), intent(in) :: x
    real (PREC), dimension(:), intent(out) :: fx

    real (PREC) :: a11, a12, a21, a22, b1, b2

    a11 = 2.0
    a12 = -5.0
    a21 = 1.5
    a22 = 2.0
    b1 = -2.0
    b2 = 3.0

    fx(1) = x(1) * a11 + x(2) * a12 - b1
    fx(2) = x(1) * a21 + x(2) * a22 - b2
end subroutine

pure subroutine func2 (x, fx)
    real (PREC), dimension(:), intent(in) :: x
    real (PREC), dimension(:), intent(out) :: fx

    fx(1) = 2.0 * x(1) ** 2 - 5 * log(x(2)) + 2
    fx(2) = 1.5 * x(1) * sqrt(x(2)) - 3.0
end subroutine

subroutine print_report (res)
    class (optim_result), intent(in) :: res

    integer, save :: ii = 1

    print "('#', t3, 'Example ', i0)", ii
    print "(t3, 'Exit status: ', i0)", res%status
    if (res%status /= OPTIM_STATUS_CONVERGED .and. len_trim(res%msg) > 0) then
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
