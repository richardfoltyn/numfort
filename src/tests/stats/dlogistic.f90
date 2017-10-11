program test_numfort_stats_dlogistic

    use iso_fortran_env
    
    use fcore_common, FC_status_t => status_t
    use fcore_testing
    use numfort_arrays
    use numfort_stats, only: logistic, dlogistic => dlogistic_real64, &
        cdf, pdf, ppf

    integer, parameter :: PREC = real64
    
    call test_all ()
    
    
contains


subroutine test_all ()
    type (test_suite) :: tests
    
    call tests%set_label ("numfort_stats_dlogistic unit tests")
    
    call test_cdf (tests)
    
    call tests%print ()
end subroutine


subroutine test_cdf (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc
    
    type (str) :: msg
    
    integer :: i
    real (PREC), dimension(3) :: mu, sigma
    real (PREC) :: diff
    real (PREC), dimension(11) :: x, fx, x2, z
    real (PREC), parameter :: tol = 1d-10
    
    tc => tests%add_test ("dlogistic CDF tests")
    
    mu = [ real(PREC) :: 0.0, 1.0, -1.0]
    sigma = [ real (PREC) :: 1.0, 0.01, 100.0]
    
    call linspace (z, - 1.0d1, 1.0d1)
    
    do i = 1, size(mu)
        x = z * sigma(i) + mu(i)
        fx = cdf (logistic, x, loc=mu(i), scale=sigma(i))
        x2 = ppf (logistic, fx, loc=mu(i), scale=sigma(i))
        
        diff = maxval(abs(x-x2))
        msg = "CDF/inv. CDF for mu=" // str(mu(i), "f0.3") // ", sigma=" &
            // str(sigma(i), "f0.3") // "; diff=" // str(diff, "e9.2")
        call tc%assert_true (diff < tol, msg)
    end do

end subroutine

end program
