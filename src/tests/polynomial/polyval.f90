

program test_polynomial_polyval

    use, intrinsic :: iso_fortran_env
    use numfort_arrays, only: vander
    use numfort_common
    use numfort_common_testing
    use numfort_polynomial, only: polyval, polyder
    
    use fcore_testing, only: test_suite, test_case
    use fcore_strings
    
    implicit none
    
    integer, parameter :: PREC = real64
    
    call test_all ()
    
    contains
    

subroutine test_all ()
    type (test_suite) :: tests 
    
    call tests%set_label ("polyfit unit tests")
    
    call test_polyval (tests)

    call tests%print ()

end subroutine


subroutine test_polyval (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc
    
    real (PREC), dimension(:), allocatable :: coefs, xp
    real (PREC) :: x, y, y_exp
    
    tc => tests%add_test ("POLYVAL unit tests")
    
    ! Degenerate "polynomial" with empty set of coefficients
    allocate (coefs(0))
    x = 0.0_PREC
    call polyval (coefs, x, y)
    call tc%assert_true (y==0.0_PREC, "Degenerate polynomial with empty coefs set")
    deallocate (coefs)

    ! Constant function, ie. polynomial of degree 0
    allocate (coefs(1), source=1.234_PREC)
    x = 0.0_PREC
    y_exp = coefs(1)
    call polyval (coefs, x, y)
    call tc%assert_true (y==y_exp, "Degree 0 polynomial")
    deallocate (coefs)

    ! Linear funcdtion,  ie. polynomial of degree 1
    allocate (coefs(2))
    coefs(1) = 1.123_PREC
    coefs(2) = 3.456_PREC
    x = 6.789_PREC
    y_exp = coefs(1) + coefs(2) * x
    call polyval (coefs, x, y)
    call tc%assert_true (all_close (y, y_exp), "Degree 1 polynomial")
    deallocate (coefs)

    ! Quadratic function
    allocate (coefs(3), source=[0.234d0, 1.234d0, 8.234d0])
    x = 0.123_PREC
    allocate (xp(3))
    call vander (x, xp, increasing=.true.)
    y_exp = dot_product (coefs, xp)
    call polyval (coefs, x, y)
    call tc%assert_true (all_close (y, y_exp), "Degree 2 polynomial")
    deallocate (coefs, xp)

    ! Cubic function
    allocate (coefs(4), source=[5.234d0, -0.234d0, 0.0d0, 1.234d0])
    allocate (xp(4))
    x = 4.234_PREC
    call vander (x, xp, increasing=.true.)
    y_exp = dot_product (coefs, xp)
    call polyval (coefs, x, y)
    call tc%assert_true (all_close (y, y_exp), "Degree 3 polynomial")
    deallocate (coefs, xp)

    ! Degree 4 polynomial
    allocate (coefs(5), source=[5.234d0, -9.234d0, 0.234d0, 0.0d0, -1.234d0])
    allocate (xp(5))
    x = -4.234_PREC
    call vander (x, xp, increasing=.true.)
    y_exp = dot_product (coefs, xp)
    call polyval (coefs, x, y)
    call tc%assert_true (all_close (y, y_exp), "Degree 4 polynomial")
    deallocate (coefs, xp)

end subroutine


end program
