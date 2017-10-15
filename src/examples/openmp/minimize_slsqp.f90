program rosenbrock_slsqp_openmp
    !*  Example program for finding minima of (parametrized) Rosenbrock 
    !   functions in parallel using OpenMP.

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env
    use numfort_arrays, only: linspace
    use numfort_optimize, optim_result => optim_result_real64, &
        workspace => workspace_real64
    
    implicit none
    
    integer, parameter :: PREC = real64
    
    call example1 ()

contains

subroutine example1 ()

    real (PREC), dimension(:,:), allocatable :: xlb, xub
    real (PREC), dimension(:,:), allocatable :: xmin, xmin_p
    real (PREC), dimension(:), allocatable :: alpha, beta, fxmin, fxmin_p
    
    real (PREC) :: NAN
    
    real (PREC) :: diff_x, diff_fx
    integer, parameter :: n = 10000
        !*  Number of problems to solve
    integer, parameter :: m = n/2
    
    NAN = ieee_value (0.0_PREC, IEEE_QUIET_NAN)
        
    allocate (alpha(n), beta(n))
    allocate (xlb(2, n), xub(2, n))
    allocate (xmin(2,n), xmin_p(2,n), fxmin(n), fxmin_p(n))
    
    call linspace (alpha, 0.8d0, 1.2d0)
    call linspace (beta, 1.2d0, 0.8d0)
    
    call linspace (xlb(1, 1:m), -1.0d0, -0.75d0)
    call linspace (xlb(2, 1:m), -0.75d0, -0.1d0)
    ! Set other of problems to have no lower bound
    xlb(:,m+1:) = NAN
    ! Set first half of problem to have no upper bound
    call linspace (xub(1, 1:m), 0.75d0, 1.0d0)
    call linspace (xub(2, 1:m), 1.0d0, 0.75d0)
    xub(:,m+1:) = NAN
    
    call solve_serial (alpha, beta, xlb, xub, xmin, fxmin)
    call solve_parallel (alpha, beta, xlb, xub, xmin_p, fxmin_p)
   
    diff_x = maxval(abs(xmin-xmin_p))
    diff_fx = maxval(abs(fxmin-fxmin_p))
    
    print '("Max. diff in xmin: ", e10.2, "; max. diff in f(xmin): ", e10.2)', &
        diff_x, diff_fx
end subroutine


subroutine solve_serial (alpha, beta, xlb, xub, xmin, fxmin)
    real (PREC), intent(in), dimension(:) :: alpha, beta
    real (PREC), intent(in), dimension(:,:) :: xlb, xub
    real (PREC), intent(out), dimension(:,:) :: xmin
    real (PREC), intent(out), dimension(:) :: fxmin

    real (PREC), dimension(2) :: args, x0
    integer :: i, n
    integer, parameter :: m = 1
    real (PREC), parameter :: tol = 1d-8
    type (optim_result) :: res
    type (workspace) :: work
    
    n = size(alpha)
    
    do i = 1, n
        args(1) = alpha(i)
        args(2) = beta(i)
        
        x0 = [1.0d-1, 1.0d-1]
        call minimize_slsqp (fobj, x0, xlb(:,i), xub(:,i), m, f_ieqcons=fconstr, &
            res=res, work=work, tol=tol, args=args)
        
        xmin(:,i) = res%x(1:2)
        fxmin(i) = res%fx(1)
    end do
end subroutine



subroutine solve_parallel (alpha, beta, xlb, xub, xmin, fxmin)
    real (PREC), intent(in), dimension(:) :: alpha, beta
    real (PREC), intent(in), dimension(:,:) :: xlb, xub
    real (PREC), intent(out), dimension(:,:) :: xmin
    real (PREC), intent(out), dimension(:) :: fxmin

    real (PREC), dimension(2) :: args, x0
    integer :: i, n
    integer, parameter :: m = 1
    real (PREC), parameter :: tol = 1d-8
    type (optim_result) :: res
    type (workspace) :: work
    
    n = size(alpha)
    
    !$omp parallel default(none) &
    !$omp private(x0, args, i, res, work) &
    !$omp shared(alpha, beta, xlb, xub, xmin, fxmin) &
    !$omp shared(n)
    
    !$omp do schedule (auto)
    do i = 1, n
        args(1) = alpha(i)
        args(2) = beta(i)
        
        x0 = [1.0d-1, 1.0d-1]
        call minimize_slsqp (fobj, x0, xlb(:,i), xub(:,i), m, f_ieqcons=fconstr, &
            res=res, work=work, tol=tol, args=args)
        
        xmin(:,i) = res%x(1:2)
        fxmin(i) = res%fx(1)
    end do
    !$omp end do
    
    call workspace_finalize (work)
    !$omp end parallel
    
end subroutine



pure subroutine fobj (x, fx, fpx, args)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), contiguous, optional :: fpx
    real (PREC), intent(in), dimension(:), optional :: args
    
    real (PREC) :: alpha, beta
    
    alpha = 1.0
    beta = 1.0
    if (present(args)) then
        alpha = args(1)
        beta = args(2)
    end if
    
    if (present(fx)) then
        ! Compute objective 
        fx = 1.0d2 * (alpha*x(2)-beta*x(1)**2.0d0)**2.0d0 + (1.0d0-x(1))**2.0d0
    end if
    
    if (present(fpx)) then
        fpx(1) = - 4.0d2 * (alpha*x(2)-beta*x(1)**2.0d0) * beta * x(1) - 2.0d0*(1.0d0-x(1))
        fpx(2) = 2.0d2 * (alpha*x(2)-x(1)**2.0d0) * alpha
    end if
end subroutine



pure subroutine fconstr (x, fx, fpx, args)
    !*  Function evaluating inequality constraints
    !   Constraints needs to be formulated such that C(x) >= 0
    real (real64), intent(in), dimension(:), contiguous :: x
    real (real64), intent(out), dimension(:), contiguous, optional :: fx
    real (real64), intent(out), dimension(:,:), contiguous, optional :: fpx
    real (real64), intent(in), dimension(:), optional :: args
    
    real (PREC) :: alpha, beta
    
    alpha = 1.0
    beta = 1.0
    if (present(args)) then
        alpha = args(1)
        beta = args(2)
    end if
 
    if (present(fx)) then
        fx = - alpha * x(1)**2.0d0 - beta * x(2) ** 2.0d0 + 1.0d0
    end if
    
    ! Jacobian of constraint function
    if (present(fpx)) then
        fpx(1,1) = -2.0d0 * x(1) * alpha
        fpx(1,2) = -2.0d0 * x(2) * beta
    end if
end subroutine



end program
