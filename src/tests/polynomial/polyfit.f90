

program test_polynomial_polyfit

    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_polynomial, only: polyfit, polyfit_deriv
    
    use fcore_testing, only: test_suite, test_case
    use fcore_strings
    
    implicit none
    
    integer, parameter :: PREC = real64
    
    call test_all ()
    
    contains
    

subroutine test_all ()
    type (test_suite) :: tests 
    
    call tests%set_label ("polyfit unit tests")
    
    call test_polyfit_deriv (tests)
    
    call tests%print ()

end subroutine


subroutine test_polyfit_deriv (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc
    
    real (PREC), dimension(:), allocatable :: y, coefs, coefs2, xp
    real (PREC) :: x
    integer :: deg, n
    
    type (status_t) :: status
    
    tc => tests%add_test ("polyfit_deriv")
    
    ! Quadratic function
    deg = 2
    n = deg + 1
    allocate (y(n), coefs(n), coefs2(n), xp(n))
    
    coefs(:) = [1.0d0, 2.0d0, 3.0d0]
    deallocate (y, coefs)
    
end subroutine


end program
