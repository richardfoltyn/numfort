
program stats_pcr_demo
    !*  Demo program for principal component analysis (PCA) and principal
    !   component regression.

    use, intrinsic :: iso_fortran_env
    use numfort
    use numfort_common
    use numfort_stats, only: std, pcr, pca, ols
    use numfort_linalg, only: inv

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()

contains

subroutine example1 ()
    !*  Compare PCR regression to OLS, using 1,..,3 principal components

    integer, parameter :: nobs = 50
    integer, parameter :: nvars = 3

    real (PREC), dimension(:,:), allocatable :: x, z, scores, loadings
    real (PREC), dimension(:), allocatable :: sval, mean_x, std_x, propvar, shocks
    real (PREC), dimension(:), allocatable :: y
    real (PREC), dimension(:, :), allocatable :: resid
        !!  Prediction errors
    real (PREC), dimension(nvars-1, nvars-1) :: rotm
    real (PREC), dimension(nvars+1, nvars+1) :: coefs
        !!  PCR coefficients for 1,2 or 3 PCR components
    real (PREC) :: theta, mean_y, std_y, std_s
    logical :: center, scale, trans_x
    integer :: i, ncomp
    type (status_t) :: status

    ! Arguments for RANDOM_SEED
    integer :: rnd_size
    integer, allocatable, dimension(:) :: rnd_put

    ! Arguments to DGEMM
    character (1) :: transa, transb
    integer :: lda, ldb, ldc, m, n, k
    real (PREC) :: alpha = 1.0_PREC, beta = 0.0_PREC

    ! Set PRNG seed
    call random_seed (size=rnd_size)
    allocate (rnd_put(rnd_size))
    rnd_put = 1

    allocate (z(nobs, nvars-1), x(nobs, nvars), y(nobs))
    call random_number (z)

    ! rescale variables to something
    z(:,1) = 1000 * z(:,1)

    theta = pi / 4.0d0
    rotm = reshape([cos(theta), sin(theta), -sin(theta), cos(theta)], &
        shape=[nvars-1, nvars-1])

    ! Rotate to get some dependence
    transa = 'N'
    transb = 'T'
    lda = nobs
    ldb = nvars - 1
    ldc = nobs
    m = nobs
    n = nvars - 1
    k = nvars - 1
    call DGEMM (transa, transb, m, n, k, alpha, z, lda, rotm, ldb, beta, x(:,1:2), ldc)

    ! Add some uncorrelated component
    call random_number (x(:, 3))
    x(:,3) = x(:,3) * 10 - 5

    ! y is a function of x w/o intercept
    call func1 (x, y)

    ! Add some noise to dependent variable
    allocate (shocks(nobs))
    call random_number (shocks)

    call std (y, m=mean_y, s=std_y, dof=1)
    ! Rescale such that std. dev. of error term is a 100th of y
    std_s = sqrt(std_y ** 2 / 1.0d4 * 12)
    shocks = (shocks - 0.5d0) * std_s
    y = y + shocks

    ! Run regular OLS
    call ols (y, x, coefs(:,1), add_const=.true., status=status)

    ! Allocate array to store residuals for OLS and for PCR with 1,2,3 components
    allocate (resid(nobs, nvars+1))
    call residuals (y, x, coefs(:,1), resid(:, 1), add_const=.true.)

    ! PCA
    center = .true.
    scale = .true.
    trans_x = .false.
    ncomp = nvars
    allocate (scores(nobs, ncomp), loadings(nvars, ncomp))
    allocate (mean_x(nvars), std_x(nvars))
    allocate (propvar(ncomp), sval(ncomp))

    call pca (x, scores, ncomp, center, scale, trans_x, sval, &
        loadings, mean_x, std_x, propvar, status)

    ! PCR: run for i = 1, 2, 3 principal components
    do i = 1, nvars
        ! First column of coefs array contains OLS estimates!
        call pcr (y, scores(:, 1:i), sval(1:i), loadings(:,1:i), coefs(:,i+1), &
            mean_x, std_x, status=status)
        ! Store in-sample prediction errors; first column is already used for OLS!
        call residuals (y, x, coefs(:, i+1), resid(:, i+1), add_const=.true.)
    end do

    call print_report (y, coefs, resid)


end subroutine

subroutine print_report (y, coefs, resid)

    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(in), dimension(:,:) :: coefs, resid

    real (PREC), dimension(size(coefs,2)) :: tmp
    real (PREC), dimension(size(coefs,2)) :: rsq, rmse

    integer :: nvars, nobs, i
    real (PREC) :: std_y, diff

    ! nvars excluding intercept
    nvars = size(coefs, 1) - 1
    nobs = size(y)

    call std (y, s=std_y, dof=1)

    print "(/, tr1, a)", "Estimated coefficients:"
    print "(/, tr5, *(a10, :, tr2))", "OLS", "PC1", "PC2", "PC3"
    do i = 1, nvars+1
        ! Store in temporary array, otherwise ifort produces runtime warnings
        tmp(1:nvars+1) = coefs(i, :)
        print "(tr5, *(e10.3, :, tr2))", tmp(1:nvars+1)
    end do

    ! OLS and PCR should yield more or less identical results if all
    ! PCs are used.
    diff = maxval(abs(coefs(:,1) - coefs(:, nvars+1)))
    print '(/, tr5, a, ": ", e10.3)', &
        "Max. abs. diff. between OLS and PCR using all comp.", diff

    print '(/, tr1, a)', "Diagnostic stats:"
    print "(/, t20, *(a10, :, tr2))", "OLS", "PC1", "PC2", "PC3"
    ! R^2
    rsq = 1 - sum(resid ** 2, dim=1) / (nobs-1) / (std_y ** 2)
    print "(tr5, a, t20, *(f10.8, :, tr2))", "R-squared", rsq
    ! RMSE = sqrt(RSS/(n-k))
    rmse = sqrt(sum(resid ** 2, dim=1) / (nobs - nvars - 1))
    print "(tr5, a, t20, *(f10.8, :, tr2))", "RMSE", rmse

    ! Print table or in-sample prediction errors for OLS and each number of PCs
    print "(/, tr1, a)", "In-sample prediction errors:"
    print "(/, tr5, *(a10, :, tr2))", "y=f(x)", "OLS", "PC1", "PC2", "PC3"
    do i = 1, nobs
        tmp = resid(i,:)
        print "(tr5, *(f10.5, :, tr2))", y(i), tmp
    end do

end subroutine

subroutine residuals (y, x, coefs, resid, add_const)

    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(in), dimension(:,:), target :: x
    real (PREC), intent(in), dimension(:), target :: coefs
    real (PREC), intent(out), dimension(:), target :: resid
    logical, intent(in), optional :: add_const

    logical :: ladd_const
    real (PREC), dimension(:,:), pointer :: ptr_x
    integer :: nobs, ncoefs, nconst

    integer :: m, n, k, lda, ldb, ldc
    real (PREC), parameter :: alpha = 1.0_PREC, beta = 0.0_PREC
    real (PREC), dimension(:,:), pointer :: ptr_coefs, ptr_resid
    character (1) :: transa, transb

    ladd_const = .true.
    if (present(add_const)) ladd_const = add_const

    nconst = 0
    if (ladd_const) nconst = nconst + 1
    nobs = size(x, 1)
    ncoefs = size(coefs)

    if (ladd_const) then
        allocate (ptr_x(nobs, ncoefs))
        ptr_x(:, 2:ncoefs) = x
        ptr_x(:,1) = 1.0_PREC
    else
        ptr_x => x
    end if

    ! Set up call to GEMM
    m = nobs
    n = 1
    k = ncoefs
    lda = nobs
    ldb = ncoefs
    ldc = nobs
    transa = 'N'
    transb = 'N'
    ptr_coefs(1:ncoefs,1:1) => coefs
    ptr_resid(1:nobs,1:1) => resid

    call DGEMM (transa, transb, m, n, k, alpha, ptr_x, lda, ptr_coefs, ldb, &
        beta, ptr_resid, ldc)

    ! Compute residuals from predicted values
    resid = y - resid

end subroutine

pure subroutine func1 (x, y)
    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(out), dimension(:) :: y

    y = sum(x**2, dim=2) / 1.0d4 + 10
end subroutine

end
