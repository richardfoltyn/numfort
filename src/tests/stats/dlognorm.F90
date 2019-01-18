

program test_numfort_stats_dlognorm

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays
    use numfort_stats, dlognorm => dlognorm_real64
    use numfort_common_testing

    use fcore_testing
    use fcore_strings

    implicit none

    integer, parameter :: PREC = real64

    type (dlognorm) :: lognorm
        ! Lognormal distribution object with default parameters

    ! Include test data generated from Scipy functions
#include "lognorm_data.F90"

    call test_all ()


    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ('Unit tests for log-normal distribution')

    call test_pdf (tests)
    call test_cdf (tests)
    call test_ppf (tests)
    call test_cdf_ppf_inv (tests)

    call tests%print ()

end subroutine



subroutine test_pdf (tests)
    !*  Unit tests for the lognormal PDF routine
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (str) :: msg

    real (PREC), dimension(:), allocatable :: fx
    real (PREC) :: sigma, mu, xi, fxi
    integer :: i, j, n
    real (PREC), parameter :: ATOL = 1.0e-10_PREC, RTOL=1.0e-10_PREC
    logical :: is_ok

    tc => tests%add_test ('Log-normal PDF tests')

    ! Compute PDF for all sigma/mu combinations used to generate the Scipy
    ! diagnostic data
    ! Number of data points for each location/scale parameter combination
    n = size(XDATA)
    allocate (fx(n))

    do i = 1, size(SCALE_PARAMS)
        sigma = SCALE_PARAMS(i)

        do j = 1, size(LOC_PARAMS)
            mu = LOC_PARAMS(j)

            ! Compare to PDF computed in Scipy
            fx(:) = pdf (lognorm, XDATA, loc=mu, scale=sigma)
            msg = 'PDF for loc =' // str(mu, 'f7.3') // '; scale =' &
                // str(sigma, 'f7.3')
            is_ok = all_close (fx, SCIPY_PDF_DATA(:,j,i), atol=ATOL, rtol=RTOL)
            call tc%assert_true (is_ok, msg)
        end do
    end do

    ! Test PDF at selected known (boundary) values
    do i = 1, size(SCALE_PARAMS)
        sigma = SCALE_PARAMS(i)

        do j = 1, size(LOC_PARAMS)
            mu = LOC_PARAMS(j)

            ! Evaluate at lower bound
            xi = 0.0
            fxi = pdf (lognorm, xi, loc=mu, scale=sigma)
            msg = 'PDF(0.0) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (fxi == 0.0, msg)

            ! Evaluate at "upper" bound
            xi = ieee_value (xi, IEEE_POSITIVE_INF)
            fxi = pdf (lognorm, xi, loc=mu, scale=sigma)
            msg = 'PDF(inf) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (fxi == 0.0, msg)

            ! Evaluate at negative values
            xi = -1.0
            fxi = pdf (lognorm, xi, loc=mu, scale=sigma)
            msg = 'PDF(-1.0) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (fxi == 0.0, msg)
        end do
    end do

end subroutine



subroutine test_cdf (tests)
    !*  Unit tests for the lognormal CDF routine
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (str) :: msg

    real (PREC), dimension(:), allocatable :: fx
    real (PREC) :: sigma, mu, xi, fxi
    integer :: i, j, n
    real (PREC), parameter :: ATOL = 1.0e-10_PREC, RTOL=1.0e-10_PREC
    logical :: is_ok

    tc => tests%add_test ('Log-normal CDF tests')

    ! Compute CDF for all sigma/mu combinations used to generate the Scipy
    ! diagnostic data
    ! Number of data points for each location/scale parameter combination
    n = size(XDATA)
    allocate (fx(n))

    do i = 1, size(SCALE_PARAMS)
        sigma = SCALE_PARAMS(i)

        do j = 1, size(LOC_PARAMS)
            mu = LOC_PARAMS(j)

            ! Compare to CDF computed in scipy
            fx(:) = cdf (lognorm, XDATA, loc=mu, scale=sigma)
            msg = 'CDF for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            is_ok = all_close (fx, SCIPY_CDF_DATA(:,j,i), atol=ATOL, rtol=RTOL)
            call tc%assert_true (is_ok, msg)
        end do
    end do

    ! Test CDF at selected known (boundary) values
    do i = 1, size(SCALE_PARAMS)
        sigma = SCALE_PARAMS(i)

        do j = 1, size(LOC_PARAMS)
            mu = LOC_PARAMS(j)

            ! Evaluate at lower bound
            xi = 0.0
            fxi = cdf (lognorm, xi, loc=mu, scale=sigma)
            msg = 'CDF(0.0) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (fxi == 0.0, msg)

            ! Evaluate at "upper" bound
            xi = ieee_value (xi, IEEE_POSITIVE_INF)
            fxi = cdf (lognorm, xi, loc=mu, scale=sigma)
            msg = 'CDF(inf) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (fxi == 1.0_PREC, msg)

            ! Evaluate at negative values
            xi = -1.0
            fxi = cdf (lognorm, xi, loc=mu, scale=sigma)
            msg = 'CDF(-1.0) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (fxi == 0.0, msg)
        end do
    end do

end subroutine



subroutine test_ppf (tests)
    !*  Unit tests for the lognormal PPF (inverse CDF) routine
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (str) :: msg

    real (PREC), dimension(:), allocatable :: xx
    real (PREC) :: sigma, mu, fxi, qi
    integer :: i, j
    real (PREC), parameter :: ATOL = 1.0e-10_PREC, RTOL=1.0e-10_PREC
    logical :: is_ok

    tc => tests%add_test ('Log-normal PPF tests')

    ! Compute PPF for all sigma/mu combinations used to generate the Scipy
    ! diagnostic data

    ! Number of data points for each location/scale parameter combination
    allocate (xx(size(QUANTILE_DATA)))

    do i = 1, size(SCALE_PARAMS)
        sigma = SCALE_PARAMS(i)

        do j = 1, size(LOC_PARAMS)
            mu = LOC_PARAMS(j)

            ! Comparte to quantiles computed in scipy
            xx(:) = ppf (lognorm, QUANTILE_DATA, loc=mu, scale=sigma)
            msg = 'PPF for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            is_ok = all_close (xx, SCIPY_PPF_DATA(:,j,i), atol=ATOL, rtol=RTOL)
            call tc%assert_true (is_ok, msg)

        end do
    end do

    ! Test PPF at selected known (boundary) values
    do i = 1, size(SCALE_PARAMS)
        sigma = SCALE_PARAMS(i)

        do j = 1, size(LOC_PARAMS)
            mu = LOC_PARAMS(j)

            ! Evaluate at lower bound
            fxi = 0.0
            qi = ppf (lognorm, fxi, loc=mu, scale=sigma)
            msg = 'PPF(0.0) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (qi == 0.0, msg)

            ! Evaluate at "upper" bound
            fxi = 1.0
            qi = ppf (lognorm, fxi, loc=mu, scale=sigma)
            msg = 'PPF(1.0) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (ieee_class (qi) == IEEE_POSITIVE_INF, msg)

            ! PPF for invalid domain
            fxi = -1.0
            qi = ppf (lognorm, fxi, loc=mu, scale=sigma)
            msg = 'PPF(-1.0) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (ieee_is_nan (qi), msg)

            fxi = 1.1
            qi = ppf (lognorm, fxi, loc=mu, scale=sigma)
            msg = 'PPF(1.1) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (ieee_is_nan (qi), msg)
        end do
    end do

end subroutine



subroutine test_cdf_ppf_inv (tests)
    ! Unit tests checking that the CDF and PPF are (approximate) inverses
    ! of each other.
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (str) :: msg

    real (PREC), parameter :: MU_ALL(*) = [-10.234d0, 0.0d0, 1.234d0]
    real (PREC), parameter :: SIGMA_ALL(*) = [0.0213d0, 1.0d0, 2.2345d0]
    real (PREC), dimension(:), allocatable :: xx, qq, xx1, qq1
    real (PREC) :: sigma, mu, diff
    integer :: i, j, n
    real (PREC), parameter :: ATOL = 1.0e-8_PREC, RTOL=1.0e-8_PREC
    logical :: is_ok

    tc => tests%add_test ('Log-normal PPF tests')

    ! Compute PPF for all sigma/mu combinations used to generate the Scipy
    ! diagnostic data

    ! Number of data points for each location/scale parameter combination
    n = 20
    allocate (xx(n), xx1(n), qq(n), qq1(n))

    call linspace (xx, 0.0_PREC, 10.0_PREC)
    call linspace (qq, 0.0_PREC, 1.0 - 1.0e-5_PREC)

    do i = 1, size(SIGMA_ALL)
        sigma = SIGMA_ALL(i)

        do j = 1, size(MU_ALL)
            mu = MU_ALL(j)

            ! Direction 1: PPF(CDF(x)) == x
            qq1(:) = cdf (lognorm, xx, loc=mu, scale=sigma)
            xx1(:) = ppf (lognorm, qq1, loc=mu, scale=sigma)
            ! Ignore points that get mapped to 1.0 as the inverse is always
            ! inf.
            where (qq1 == 0.0_PREC .or. qq1 == 1.0_PREC)
                xx1 = xx
            end where

            diff = maxval(abs(xx1-xx))

            is_ok = all_close (xx1, xx, atol=ATOL, rtol=RTOL)
            msg = 'PPF(CDF(x)) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (is_ok, msg)

            ! Direction 2: CDF(PPF(q)) == q
            xx1(:) = ppf (lognorm, qq, loc=mu, scale=sigma)
            qq1(:) = cdf (lognorm, xx1, loc=mu, scale=sigma)
            is_ok = all_close (qq1, qq, atol=ATOL, rtol=RTOL)
            msg = 'CDF(PPF(q)) for loc =' // str(mu, 'f7.3') // '; scale = ' &
                // str(sigma, 'f7.3')
            call tc%assert_true (is_ok, msg)

        end do
    end do

end subroutine

end program
