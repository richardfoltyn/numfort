program rosenbrock_slsqp_openmp
    !*  Example program for finding minima of (parametrized) Rosenbrock
    !   functions in parallel using OpenMP.

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env
    use omp_lib
    use numfort_arrays, only: linspace
    use numfort_optimize, optim_result => optim_result_real64, &
        workspace => workspace_real64

    implicit none

    integer, parameter :: PREC = real64

    type, extends(args_data) :: args_t
       real (PREC) :: alpha, beta
    end type

    interface dynamic_cast
        procedure cast_to_args
    end interface

    call example1 ()

contains

subroutine example1 ()

    real (PREC), dimension(:,:), allocatable :: xlb, xub
    real (PREC), dimension(:,:), allocatable :: xmin, xmin_p
    real (PREC), dimension(:), allocatable :: alpha, beta, fxmin, fxmin_p

    real (PREC) :: NAN, POSITIVE_INF, NEGATIVE_INF

    real (PREC) :: diff_x, diff_fx
    integer, parameter :: n = 100000
        !*  Number of problems to solve
    integer, parameter :: m = n/2

    NAN = ieee_value (0.0_PREC, IEEE_QUIET_NAN)
    POSITIVE_INF = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    NEGATIVE_INF = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)

    allocate (alpha(n), beta(n))
    allocate (xlb(2, n), xub(2, n))
    allocate (xmin(2,n), xmin_p(2,n), fxmin(n), fxmin_p(n))

    call linspace (alpha, 0.8d0, 1.2d0)
    call linspace (beta, 1.2d0, 0.8d0)

    call linspace (xlb(1, 1:m), -1.0d0, -0.75d0)
    call linspace (xlb(2, 1:m), -0.75d0, -0.1d0)
    ! Set other of problems to have no lower bound
    xlb(:,m+1:) = NEGATIVE_INF
    ! Set first half of problem to have no upper bound
    call linspace (xub(1, 1:m), 0.75d0, 1.0d0)
    call linspace (xub(2, 1:m), 1.0d0, 0.75d0)
    xub(:,m+1:) = POSITIVE_INF

    call solve_serial (alpha, beta, xlb, xub, xmin, fxmin)
    call solve_parallel (alpha, beta, xlb, xub, xmin_p, fxmin_p)

    diff_x = maxval(abs(xmin-xmin_p))
    diff_fx = maxval(abs(fxmin-fxmin_p))

    print *, "Comparing serial to parallel results"
    print '(tr3, "Max. diff in xmin: ", e10.2, "; max. diff in f(xmin): ", e10.2)', &
        diff_x, diff_fx
end subroutine


subroutine solve_serial (alpha, beta, xlb, xub, xmin, fxmin)
    real (PREC), intent(in), dimension(:) :: alpha, beta
    real (PREC), intent(in), dimension(:,:) :: xlb, xub
    real (PREC), intent(out), dimension(:,:) :: xmin
    real (PREC), intent(out), dimension(:) :: fxmin

    real (PREC), dimension(2) :: x0
    integer :: i, n
    integer, parameter :: m = 1
    real (PREC), parameter :: tol = 1d-8
    type (optim_result) :: res
    type (workspace) :: work
    type (args_t) :: args

    real :: tic, toc

    n = size(alpha)

    call cpu_time (tic)

    do i = 1, n
        args%alpha = alpha(i)
        args%beta = beta(i)

        x0 = [1.0d-1, 1.0d-1]
        call minimize_slsqp (fobj, x0, args, xlb(:,i), xub(:,i), m, f_ieqcons=fconstr, &
            res=res, work=work, tol=tol)

        xmin(:,i) = res%x(1:2)
        fxmin(i) = res%fx(1)
    end do

    call cpu_time (toc)

    print *, "Serial execution"
    print '(tr3, "Time elapsed: ", f0.2, " seconds")', toc-tic
end subroutine



subroutine solve_parallel (alpha, beta, xlb, xub, xmin, fxmin)
    real (PREC), intent(in), dimension(:) :: alpha, beta
    real (PREC), intent(in), dimension(:,:) :: xlb, xub
    real (PREC), intent(out), dimension(:,:) :: xmin
    real (PREC), intent(out), dimension(:) :: fxmin

    real (PREC), dimension(2) :: x0
    integer :: i, n
    integer, parameter :: m = 1
    real (PREC), parameter :: tol = 1d-8
    type (optim_result) :: res
    type (workspace) :: work
    integer :: ncpus
    type (args_t) :: args

    real (PREC) :: tic, toc

    n = size(alpha)

    ! Use OMP timing functions instead of CPU_TIME, as the latter returns
    ! CPU time summed across all CPUs.
    tic = omp_get_wtime ()

    !$omp parallel default(none) &
    !$omp private(x0, args, i, res, work) &
    !$omp shared(alpha, beta, xlb, xub, xmin, fxmin) &
    !$omp shared(n, ncpus)

    ncpus = omp_get_num_threads ()

    !$omp do schedule (auto)
    do i = 1, n
        args%alpha = alpha(i)
        args%beta = beta(i)

        x0 = [1.0d-1, 1.0d-1]
        call minimize_slsqp (fobj, x0, args, xlb(:,i), xub(:,i), m, f_ieqcons=fconstr, &
            res=res, work=work, tol=tol)

        xmin(:,i) = res%x(1:2)
        fxmin(i) = res%fx(1)
    end do
    !$omp end do

    !$omp end parallel

    toc = omp_get_wtime ()

    print *, "Parallel execution"
    print '(tr3, "Number of cores: ", i0)', ncpus
    print '(tr3, "Time elapsed: ", f0.2, " seconds")', toc-tic
end subroutine



subroutine fobj (x, args, fx, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), contiguous, optional :: fpx

    real (PREC) :: alpha, beta
    type (args_t), pointer :: largs

    call dynamic_cast (args, largs)

    alpha = largs%alpha
    beta = largs%beta

    if (present(fx)) then
        ! Compute objective
        fx = 1.0d2 * (alpha*x(2)-beta*x(1)**2.0d0)**2.0d0 + (1.0d0-x(1))**2.0d0
    end if

    if (present(fpx)) then
        fpx(1) = - 4.0d2 * (alpha*x(2)-beta*x(1)**2.0d0) * beta * x(1) - 2.0d0*(1.0d0-x(1))
        fpx(2) = 2.0d2 * (alpha*x(2)-x(1)**2.0d0) * alpha
    end if
end subroutine



subroutine fconstr (x, args, fx, fpx)
    !*  Function evaluating inequality constraints
    !   Constraints needs to be formulated such that C(x) >= 0
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx

    real (PREC) :: alpha, beta
    type (args_t), pointer :: largs

    call dynamic_cast (args, largs)

    alpha = largs%alpha
    beta = largs%beta

    if (present(fx)) then
        fx = - alpha * x(1)**2.0d0 - beta * x(2) ** 2.0d0 + 1.0d0
    end if

    ! Jacobian of constraint function
    if (present(fpx)) then
        fpx(1,1) = -2.0d0 * x(1) * alpha
        fpx(1,2) = -2.0d0 * x(2) * beta
    end if
end subroutine


subroutine cast_to_args (tgt, ptr)
    class (args_data), intent(in), target :: tgt
    type (args_t), intent(inout), pointer :: ptr

    nullify (ptr)

    select type (tgt)
    type is (args_t)
        ptr => tgt
    end select
end subroutine


end program
