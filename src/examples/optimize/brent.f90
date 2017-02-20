program brent

    use numfort_optimize
    use iso_fortran_env

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()

contains


subroutine example1 ()

    real (PREC) :: a, b, x0
    type (optim_result) :: res

    a = -10
    b = (-3-sqrt(5.0d0))/2.0 + 0.1

    call brentq (func1, a, b, x0=x0, res=res)
    call print_report (res)

    a = b
    b = 1.9d0

    call brentq (func1, a, b, x0=x0, res=res)
    call print_report (res)

    a = x0 + 0.01
    b = 100.0

    call brentq (func1, a, b, x0=x0, res=res)
    call print_report (res)

    ! test with invalid bracketing interval
    a = 2.1d0

    call brentq (func1, a, b, x0=x0, res=res)
    call print_report (res)

end subroutine

function func1 (x) result (fx)
    real (PREC), intent(in) :: x
    real (PREC) :: fx

    ! use 3rd-degree polynomial with roots at (-3-sqrt(5))/2,
    ! (-3+sqrt(5))/2 and 2
    fx = (x+3) * (x-1)**2 - 5
end function

subroutine print_report (res)
    class (optim_result), intent(in) :: res

    integer, save :: ii = 1

    print "('#', t3, 'Example ', i0)", ii
    print "(t3, 'Exit status: ', i0)", res%status
    if ((.not. res%success) .and. len_trim(res%msg) > 0) then
        print "(t3, 'Message: ', a)", res%msg
    end if
    if (res%success) then
        print "(t3, 'Root located at: ', g23.15e2)", res%x
        print "(t3, 'Function value at root: ', g23.15e2)", res%fx
        print "(t3, 'Number of iterations: ', i0)", res%nit
        print "(t3, 'Number of function evaluations: ', i0)", res%nfev
    endif

    ii = ii + 1

end subroutine

end
