
program test_core_combinatorics

    use, intrinsic :: iso_fortran_env
    use numfort_arrays
    use numfort_common
    use numfort_common_testing
    use numfort_core

    use fcore_testing, only: test_suite, test_case
    use fcore_strings

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains


subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_core combinatorics unit tests")

    call test_factorial (tests)
    call test_comb (tests)

    call tests%print ()

end subroutine


subroutine test_factorial (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    integer, parameter :: N = 10
    integer, parameter, dimension(N) :: fact_true = &
        [ 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880]
    integer, dimension(N) :: x, fx

    tc => tests%add_test ("Factorial function")

    ! Test factorial function for first 10 values
    call arange(x, 0)
    fx = factorial (x)

    call tc%assert_true (all(fx == fact_true), &
        "Factorial function for values [0,...,9]")


end subroutine


subroutine test_comb (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (str) :: msg

    integer, PARAMETER :: NMAX = 10
    logical :: ok
    integer (int64) :: res, res_true, n, k
        ! Note: need to use in64 otherwise the results for higher N won't fit
        ! into default INTEGER!

    tc => tests%add_test ("COMB with repetition")

    do n = 0, NMAX
        ok = .true.
        do k = 0, n
            res = comb (n, k, repetition=.true.)
            res_true = factorial (n+k-1) / factorial(k) / factorial (n-1)
            ok = ok .and. (res == res_true)
        end do

        msg = "COMB with N=" // str(n) // "; k=[0,...," // str(n) // "]"
        call tc%assert_true (ok, msg)
    end do

    tc => tests%add_test ("COMB without repetition")

    do n = 0, NMAX
        ok = .true.
        do k = 0, n
            res = comb (n, k, repetition=.false.)
            res_true = factorial (n) / factorial(k) / factorial (n-k)
            ok = ok .and. (res == res_true)
        end do

        msg = "COMB with N=" // str(n) // "; k=[0,...," // str(n) // "]"
        call tc%assert_true (ok, msg)
    end do


end subroutine

end