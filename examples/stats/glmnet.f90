

program example_numfort_stats_glmnet

    use, intrinsic :: iso_fortran_env

    use numfort_stats, dnorm => dnorm_real64
    use numfort_stats_lm, lm_result => lm_result_real64, lm_config => lm_config_real64
    use numfort_stats_glmnet, enet_config => enet_config_real64, &
        enet_result => enet_result_real64

    use blas95, only: BLAS_GEMV => GEMV, BLAS_GEMM => GEMM

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()
    call example2 ()

    contains


subroutine example1 ()
    !*  EXAMPLE1 demonstrates the usage of elastic net and compares it to
    !   OLS for a regressor matrix with many correlated regressors.
    real (PREC), dimension(:,:), allocatable :: X
    real (PREC), dimension(:), allocatable :: y, coefs_true, coefs_ols
    real (PREC), dimension(:), allocatable :: rwork
    real (PREC), dimension(:), allocatable :: mean_x, std_x, beta, eps
    real (PREC) :: rmse, rsq_enet, rsq_ols
    real (PREC) :: intercept_true, intercept_ols
    real (PREC) :: alpha, l1_ratio
    integer :: Nobs, Nvars, i, k
    type (status_t) :: status
    type (enet_result) :: res
    type (enet_config) :: conf
    type (lm_result) :: lm
    type (lm_config) :: conf_lm
    character (*), parameter :: FMT_COEF = '(*(g14.3,:, " "))'

    call set_seed (1234)

    ! --- Create data ---

    Nvars = 20
    k = Nvars/2 + 1
    Nobs = 1000

    allocate (X(Nobs, Nvars))
    allocate (mean_x(Nvars), std_x(Nvars))

    call random_number (mean_x(1:k))
    call random_number (std_x(1:k))

    mean_x(1:k) = (mean_x(1:k) - 0.50d0) * 2.0
    std_x(1:k) = std_x(1:k) * 2.0d0

    do i = 1, k
        call rvs (norm, X(:,i), loc=mean_x(i), scale=std_x(i))
    end do

    ! Construct remaining variables as functions of first 1:k
    allocate (eps(Nobs))
    allocate (beta(k))
    do i = k+1, Nvars
        call rvs (norm, eps, scale=0.1_PREC)
        call random_number (beta)
        beta(:) = (beta - 0.5) * 2.0d0
        call BLAS_GEMV (X(:,1:k), beta, X(:,i))
        X(:,i) = X(:,i) + eps
    end do

    ! --- True coefs ---

    allocate (coefs_true(Nvars))
    do i = 1, Nvars
        coefs_true(i) = (-1.0_PREC)**i * exp(-i/10.0_PREC)
    end do

    ! Set half to 0
    coefs_true(Nvars/2+1:) = 0.0

    ! True intercept
    intercept_true = 1.234d0

    ! --- Outcome variable ---

    allocate (y(Nobs), source=intercept_true)
    call BLAS_GEMV (X, coefs_true, y, beta=1.0_PREC)

    ! --- OLS ---

    allocate (coefs_ols(Nvars))
    conf_lm = lm_config(add_intercept=.true., rcond=0.0_PREC)
    call ols (conf_lm, X, y, coefs_ols, intercept=intercept_ols, res=lm, status=status)
    call post_estim (lm, X, y, rsq=rsq_ols)

    ! --- CV elastic net ---
    l1_ratio = 0.9

    conf%cv_n = min(100, Nobs/10)
    conf%alpha_eps = 1.0e-5_PREC
    conf%cv_parsimonious = .true.
    call enet_cv (conf, X, y, alpha, l1_ratio, rmse=rmse, status=status)

    ! Fit coefs using CV'd alpha
    call enet_fit (conf, X, y, alpha, l1_ratio, res=res, status=status)
    call post_estim (res, X, y, rsq=rsq_enet)

    ! Print OLS and ENET coefs

    print '(*(a15,:," "))', 'True         ', 'OLS          ', 'Elastic Net  '
    print FMT_COEF, intercept_true, intercept_ols, res%intercept
    do i = 1, Nvars
        print FMT_COEF, coefs_true(i), coefs_ols(i), res%coefs(i)
    end do

    print '(*(a15,:," "))', 'R-squared    ', 'OLS          ', 'Elastic Net  '
    print '(tr15, f15.10, f15.10)', rsq_ols, rsq_enet

end subroutine



subroutine example2 ()
    !*  EXAMPLE2 demonstrates the usage of elastic net and compares it to
    !   OLS for a regressor matrix with many correlated regressors,
    !   for multiple outcome variables
    real (PREC), dimension(:,:), allocatable :: X
    real (PREC), dimension(:,:), allocatable :: Y, coefs_true, coefs_ols
    real (PREC), dimension(:,:), allocatable :: rwork
    real (PREC), dimension(:), allocatable :: mean_x, std_x, beta, eps
    real (PREC), dimension(:), allocatable :: intercept_true, intercept_ols
    real (PREC) :: rmse, rsq_enet, rsq_ols
    real (PREC) :: alpha, l1_ratio
    integer :: Nobs, Nvars, Nlhs, i, j, k
    type (status_t) :: status
    type (enet_result) :: res
    type (enet_config) :: conf
    type (lm_result) :: lm
    type (lm_config) :: conf_lm
    character (*), parameter :: FMT_COEF = '(*(f9.4,:, " "))'

    call set_seed (1234)

    ! --- Create data ---

    Nvars = 20
    k = Nvars/2 + 1
    Nobs = 1000
    Nlhs = 2

    allocate (X(Nobs, Nvars))
    allocate (mean_x(Nvars), std_x(Nvars))

    call random_number (mean_x(1:k))
    call random_number (std_x(1:k))

    mean_x(1:k) = (mean_x(1:k) - 0.50d0) * 2.0
    std_x(1:k) = std_x(1:k) * 2.0d0

    do i = 1, k
        call rvs (norm, X(:,i), loc=mean_x(i), scale=std_x(i))
    end do

    ! Construct remaining variables as functions of first 1:k
    allocate (eps(Nobs))
    allocate (beta(k))
    do i = k+1, Nvars
        call rvs (norm, eps, scale=0.1_PREC)
        call random_number (beta)
        beta(:) = (beta - 0.5) * 2.0d0
        call BLAS_GEMV (X(:,1:k), beta, X(:,i))
        X(:,i) = X(:,i) + eps
    end do

    ! --- True coefs ---

    allocate (coefs_true(Nvars,Nlhs))
    do j = 1, Nlhs
        do i = 1, Nvars
            coefs_true(i,j) = (-1.0_PREC)**j * (-1.0_PREC)**i * exp(-i/10.0_PREC) * j
        end do

        ! Set fraction to 0
        coefs_true(Nvars/j:,j) = 0.0
    end do

    ! True intercept
    allocate (intercept_true(Nlhs))
    call rvs (norm, intercept_true, scale=2.0_PREC)

    ! --- Outcome variable ---

    allocate (Y(Nobs,Nlhs))
    Y(:,:) = spread(intercept_true, dim=1, ncopies=Nobs)
    call BLAS_GEMM (X, coefs_true, Y, beta=1.0_PREC)

    ! --- OLS ---

    allocate (coefs_ols(Nvars,Nlhs), intercept_ols(Nlhs))
    conf_lm = lm_config(add_intercept=.true., rcond=0.0_PREC)

    call ols (conf_lm, X, Y, coefs_ols, intercept=intercept_ols, status=status)

    ! --- CV elastic net ---
    l1_ratio = 0.9

    conf%cv_n = min(100, Nobs/10)
    conf%alpha_eps = 1.0e-5_PREC
    conf%cv_parsimonious = .true.
    call enet_cv (conf, X, Y, alpha, l1_ratio, rmse=rmse, status=status)

    ! Fit coefs using CV'd alpha
    call enet_fit (conf, X, Y, alpha, l1_ratio, res=res, status=status)

    ! Print OLS and ENET coefs

    print '(*(a20,:," "))', 'True         ', 'OLS          ', 'Elastic Net  '
    print FMT_COEF, intercept_true, intercept_ols, res%intercept_multi
    do i = 1, Nvars
        print FMT_COEF, coefs_true(i,:), coefs_ols(i,:), res%coefs_multi(i,:)
    end do

end subroutine


end
