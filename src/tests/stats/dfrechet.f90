program test_numfort_stats_dfrechet

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env
    
    use fcore_common, FC_status_t => status_t
    use fcore_testing
    use numfort_common_testing, only: all_close
    use numfort_arrays
    use numfort_stats, only: dfrechet => dfrechet_real64, &
        cdf, pdf, ppf
        
    implicit none

    integer, parameter :: PREC = real64
    
    call test_all ()
    
contains


subroutine test_all ()
    type (test_suite) :: tests
    
    call tests%set_label ("numfort_stats_dfrechet unit tests")
    
    call test_pdf (tests)
    call test_cdf (tests)
    call test_ppf (tests)
    
    call tests%print ()
end subroutine


subroutine test_pdf (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    type (str) :: msg
    type (dfrechet) :: frechet
    real (PREC) :: x, fx, shp, desired
    logical :: res
    
    tc => tests%add_test ("dfrechet PDF tests")
    
    ! Compute to values calculated in Mathematica, taken from Scipy project
    x = 2.5d0
    shp = 0.5d0
    desired = 0.067202904517205343834d0
    fx = pdf (frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (fx, desired, rtol=1d-12)
    msg = "PPF(" // str(x, "f0.3") // "; shape=" // str(shp, "f0.3") &
        // ") == Mathematica value"
    call tc%assert_true (res, msg)
    
    x = 3.0d0
    shp = 13.0d0
    desired = 2.7179753508900831584d-6
    fx = pdf (frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (fx, desired, rtol=1d-12)
    msg = "PPF(" // str(x, "f0.3") // "; shape=" // str(shp, "f0.3") &
        // ") == Mathematica value"
    call tc%assert_true (res, msg)
end subroutine




subroutine test_cdf (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc
    
    type (str) :: msg
    
    integer, parameter :: NX = 11
    integer :: i, j, k
    real (PREC), dimension(100) :: m, s, alpha
    real (PREC) :: diff
    real (PREC), dimension(NX) :: xx, fxx, xx2, z, diff_x
    real (PREC) :: x, fx
    real (PREC), parameter :: tol = 1d-10
    real (PREC) :: desired, shp
    real (PREC), dimension(100) :: power_grid
    real (PREC), dimension(NX) :: power_grid2
    logical, dimension(NX) :: mask
    logical :: res
    
    type (dfrechet) :: frechet
    
    tc => tests%add_test ("dfrechet CDF tests")
   
    call linspace (power_grid, 0.0d0, 1.0d0)
    power_grid = power_grid ** 3
    
    call linspace (power_grid2, 0.0d0, 1.0d0)
    power_grid2 = power_grid2 ** 3
    
    m = power_grid * 1.0d5
    s = 1.0d-4 + power_grid * 1.0d3
    alpha = 1.0d-4 + power_grid * 1.0d3
    
    z = 1d-5 + power_grid2 * 1.0d2
   
    diff = 0

!     do i = 1, size(m)
!         do j = 1, size(s)
!             do k = 1, size(alpha)
!                 x = z * s(j) + m(i)
!                 fx = cdf (frechet, x, loc=m(i), scale=s(j), shape=alpha(k))
!                 x2 = ppf (frechet, fx, loc=m(i), scale=s(j), shape=alpha(k))
!                 
!                 mask = ieee_is_finite (x2)
!                 
!                 where (mask)
!                     diff_x = abs(x-x2)
!                 else where
!                     diff_x = 0
!                 end where
!                 
!                 diff = max(diff, maxval(diff_x))
!             end do
!         end do
!     end do
!     
!     msg = "CDF/inv. CDF for range of loc/scale/shape params; max diff.=" & 
!         // str(diff, 'e9.2')
!     call tc%assert_true (diff < tol, msg)

    ! Compare to values calculated in Mathematica, taken from the Scipy project
    x = 2.5d0
    shp = 0.5d0
    fx = cdf(frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    desired = 0.53128560913296781152d0
    res = all_close (fx, desired, rtol=1d-12)
    msg = "CDF(" // str(x, "f0.3") // "; shape=" // str(shp, "f0.3") &
            // ") == Mathematica value"
    call tc%assert_true (res, msg)
    
    x = 1d-3
    fx = cdf(frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    desired = 1.84672666240969310047d-14
    res = all_close (fx, desired, rtol=1d-12)
    msg = "CDF(" // str(x, "f0.3") // "; shape=" // str(shp, "f0.3") &
            // ") == Mathematica value"
    call tc%assert_true (res, msg)
    
    x = 3.0d0
    shp = 13.0d0
    fx = cdf(frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    desired = 0.99999937277472231955d0
    res = all_close (fx, desired, rtol=1d-12)
    msg = "CDF(" // str(x, "f0.3") // "; shape=" // str(shp, "f0.3") &
            // ") == Mathematica value"
    call tc%assert_true (res, msg)

end subroutine


subroutine test_ppf (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    type (str) :: msg
    type (dfrechet) :: frechet
    real (PREC) :: q, x, shp, desired
    logical :: res
    
    tc => tests%add_test ("dfrechet PPF tests")
    
    ! Compute to values calculated in Mathematica, taken from Scipy project
    q = 0.975d0
    shp = 0.5d0
    desired = 1560.0833306626037620d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    msg = "PPF(" // str(q, "f0.3") // "; shape=" // str(shp, "f0.3") &
        // ") == Mathematica value"
    call tc%assert_true (res, msg)
    
    q = 0.5d0
    shp = 0.5d0
    desired = 2.0813689810056077979d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    msg = "PPF(" // str(q, "f0.3") // "; shape=" // str(shp, "f0.3") &
        // ") == Mathematica value"
    call tc%assert_true (res, msg)
    
    q = 0.975d0
    shp = 13.0d0
    desired = 1.3268241778612848659d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    msg = "PPF(" // str(q, "f0.3") // "; shape=" // str(shp, "f0.3") &
        // ") == Mathematica value"
    call tc%assert_true (res, msg)
    
    q = 0.5d0
    shp = 13.0d0
    desired = 1.0285944941498834248d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    msg = "PPF(" // str(q, "f0.3") // "; shape=" // str(shp, "f0.3") &
        // ") == Mathematica value"
    call tc%assert_true (res, msg)
    
    q = 0.01d0
    shp = 13.0d0
    desired = 0.88916242414515439302d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    msg = "PPF(" // str(q, "f0.3") // "; shape=" // str(shp, "f0.3") &
        // ") == Mathematica value"
    call tc%assert_true (res, msg)
    
    q = 1d-14
    shp = 13.0d0
    desired =  0.76554999794716080164d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    msg = "PPF(" // str(q, "e8.2") // "; shape=" // str(shp, "f0.3") &
        // ") == Mathematica value"
    call tc%assert_true (res, msg)
end subroutine

end program
