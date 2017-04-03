! Sample code demonstrating the use of curfit and splev from
! the numfort_interpolate module, using OpenMP
! Author: Richard Foltyn

program curfit_demo_omp

    use iso_fortran_env
    use numfort_interpolate, workspace => workspace_real64
    use omp_lib

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()

contains

subroutine example1 ()

    ! number of points
    integer, parameter :: m = 25
    integer, parameter :: ns = 10

    real (PREC), dimension(m) :: x, y, eps, w
    real (PREC), dimension(m, ns) :: yhat
    real (PREC), dimension(ns) :: svec
    integer :: i, k, nest, iopt, nseed, ext
    real (PREC), dimension(:), allocatable :: seed
    real (PREC), dimension(:, :), allocatable :: knots, coefs
    integer, dimension(ns) :: nknots
    type (status_t), dimension(ns) :: status
    real (PREC), dimension(ns) :: ssr

    ! arrays to hold sequential results
    real (PREC), dimension(:, :), allocatable :: knots2, coefs2
    real (PREC), dimension(ns) :: ssr2
    integer, dimension(ns) :: nknots2
    type (status_t), dimension(ns) :: status2
    real (PREC), dimension(:,:), allocatable :: yhat2

    ! locals private to each thread
    integer :: tid, j
    real (PREC) :: s
    type (workspace), pointer :: ws

    ! create vector of smoothness parameters
    svec = 1.0d1 ** [(i, i=3, 3-ns+1, -1)]

    x = [real (PREC) :: (i, i = 0, m-1)]
    ! create some reproducible random numbers
    call random_seed (size=nseed)
    allocate (seed(nseed))
    seed = 1.0_PREC
    call random_number (eps)
    eps = (eps - 0.5d0) * 2
    ! create something that looks similar to a utility function, plus add
    ! some noise
    y = log(x + 0.1) + eps
    ! spline degree
    k = 3
    iopt = 0
    w = 1/abs(eps)
    ! evaluate spline at original x points
    ext = NF_INTERP_EVAL_BOUNDARY

    nest = curfit_get_nest (m=m, k=k)
    allocate (knots(nest, ns), coefs(nest, ns))

    !$omp parallel default(none) private(j, ws, s, tid) &
    !$omp shared(knots2, coefs2, nknots2, ssr2, status2) &
    !$omp shared(svec, x, y, k, knots, coefs, nknots, iopt, w, ssr, status, nest) &
    !$omp shared(yhat, ext)
    tid = omp_get_thread_num ()
    allocate (ws)

    !$omp do
    do j = 1, ns
        s = svec(j)
        call fit (x, y, k, s, knots(:, j), coefs(:, j), nknots(j), &
            ws, ssr(j), status(j), yhat(:, j), iopt, w, ext=ext)
    end do
    !$omp end do
    deallocate (ws)
    !$omp end parallel

    do i = 1, ns
       call print_report (iopt, svec(i), k, ssr(i), status(i), &
           nknots(i), knots(:, i), coefs(:, i), x, y, yhat(:, i), i)
    end do

    ! compute sequentially
    allocate (knots2(nest, ns), coefs2(nest, ns), yhat2(m, ns))
    allocate (ws)

    do j = 1, ns
        s = svec(j)
        call fit (x, y, k, s, knots2(:, j), coefs2(:, j), nknots2(j), &
            ws, ssr2(j), status2(j), yhat2(:, j), iopt, w, ext=ext)
    end do

    deallocate (ws)

    print *, "Summary of OpenMP vs. sequential results"
    do j = 1, ns
        print '(i3, a, es9.2e2, "; ", a, es9.2e2, "; ", a, l1, "; ", a, es9.2e2, "; ", a, l1)', j, &
            "d(knots)=", maxval(abs(knots(:,j)-knots2(:,j))), &
            "d(coefs)=", maxval(abs(coefs(:,j)-coefs2(:,j))), &
            "nknots eq.=", nknots(j) == nknots2(j), &
            "d(ssr)=", abs(ssr(j)-ssr2(j)), &
            "status eq.=", status(j) == status2(j)
    end do

end subroutine

pure subroutine fit (x, y, k, s, knots, coefs, nknots, ws, ssr, status, yhat, &
        iopt, w, ext)
    real (PREC), intent(in), dimension(:) :: x, y
    integer, intent(in) :: k
    real (PREC), intent(in) :: s
    real (PREC), intent(out), dimension(:) :: knots, coefs
    integer, intent(out) :: nknots
    type (workspace), intent(in out) :: ws
    real (PREC), intent(out) :: ssr
    type (status_t), intent(out) :: status
    real (PREC), intent(out), dimension(:) :: yhat
    integer, intent(in), optional :: iopt
    real (PREC), intent(in), dimension(:), optional :: w
    integer (NF_ENUM_KIND), intent(in), optional :: ext

    call curfit (x, y, k, s, knots, coefs, nknots, &
       iopt=iopt, work=ws, w=w, ssr=ssr, status=status)

    if (.not. (NF_STATUS_INVALID_ARG .in. status)) then
       call splev (knots, coefs, nknots, k, x, yhat, ext=ext)
    end if
end subroutine

subroutine print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, yhat, counter)
    integer, intent(in) :: iopt, k, n, counter
    type (status_t), intent(in) :: status
    real (PREC), intent(in) :: s, ssr, knots(:), coefs(:)
    real (PREC), intent(in), dimension(:) :: x, y, yhat

    integer :: i, nstatus
    integer, dimension(NF_MAX_STATUS_CODES) :: istatus

    call status_decode (status, istatus, nstatus)

    print "(/,'(', i0, ')', t6, 'iopt: ', i2, '; smoothing factor: ', es12.5e2, '; spline degree: ', i1)", &
        counter, iopt, s, k
    print "(t6, 'status code: ', *(i0, :, ', '))", istatus(1:nstatus)
    ! Approximation is returned even if status /= NF_STATUS_OK as long as
    ! input arguments were valid.
    if (.not. (NF_STATUS_INVALID_ARG .in. status)) then
        print "(t6, 'SSR: ', es12.5e2)", ssr
        print "(t6, 'Number of knots: ', i0)", n
        print "(t6, 'Knots: ', *(t14, 8(f8.3, :, ', '), :, /))", knots(1:n)
        print "(t6, 'Coefs: ', *(t14, 8(f8.3, :, ', '), :, /))", coefs(1:n)
        print "(t6, a)", 'Evaluated points:'
        print "(t6, 5(3(a5, tr1), tr2))", ('x(i)', 'y(i)', 'sp(i)', i=1,5)
        print "(*(t6, 5(3(f5.1, :, tr1), tr2), :, /))", (x(i), y(i), yhat(i), i=1,size(x))
    end if
end subroutine

end program
