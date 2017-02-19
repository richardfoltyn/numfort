! Sample code demonstrating the use of curfit and splev from
! the numfort_interpolate module, using OpenMP
! Author: Richard Foltyn

program curfit_demo_omp

    use iso_fortran_env
    use numfort_common, only: decode_status
    use numfort_interpolate
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
    real (PREC), dimension(10) :: svec
    integer :: i, k, nest, iopt, nseed, ext
    real (PREC), dimension(:, :), allocatable :: knots, coefs
    real (PREC), dimension(:), allocatable :: seed
    integer, dimension(:), allocatable :: nknots
    real (PREC), dimension(:), allocatable :: knots_i, coefs_i, yhat_i
    integer :: nknots_i
    integer (NF_ENUM_KIND), dimension(:), allocatable :: status
    integer (NF_ENUM_KIND) :: status_i
    logical :: stat_ok
    real (PREC), dimension(:), allocatable :: ssr
    real (PREC) :: ssr_i

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
    eps = eps * 2
    ! create something that looks similar to a utility function, plus add
    ! some noise
    y = log(x + 0.1) + eps
    ! spline degree
    k = 3
    iopt = 0
    w = 1/eps
    ! evaluate spline at original x points
    ext = NF_INTERP_EVAL_BOUNDARY

    nest = curfit_get_nest (m=m, k=k)
    allocate (knots(nest, ns), coefs(nest, ns), nknots(ns))
    allocate (ssr(ns), status(ns))

    !$omp parallel default(none) private(j, ws, s, stat_ok, tid) &
    !$omp private(knots_i, coefs_i, nknots_i, ssr_i, status_i, yhat_i) &
    !$omp shared(svec, x, y, k, knots, coefs, nknots, iopt, w, ssr, status, nest) &
    !$omp shared(yhat, ext)
    tid = omp_get_thread_num ()
    allocate (ws)
    allocate (knots_i(nest), coefs_i(nest), yhat_i(m))
    !$omp do
    do j = 1, ns
        s = svec(j)
        call curfit (x, y, k, s, knots(:,j), coefs(:,j), nknots(j), &
           iopt=iopt, work=ws, ssr=ssr(j), status=status(j))

        stat_ok = (status(j) == NF_STATUS_OK) .or. &
           (status(j) /= NF_STATUS_INVALID_ARG .and. ssr(j) < s)

        if (status(j) /= NF_STATUS_INVALID_ARG) then
           call splev (knots(:,j), coefs(:,j), nknots(j), k, x, yhat(:, j), ext)
        end if

        ! call curfit (x, y, k, s, knots_i, coefs_i, nknots_i, iopt=iopt, &
        ! work=ws, ssr=ssr_i, status=status_i)
        !
        ! stat_ok = (status_i == NF_STATUS_OK) .or. &
        !    (status_i /= NF_STATUS_INVALID_ARG .and. ssr_i < s)
        !
        ! if (status_i /= NF_STATUS_INVALID_ARG) then
        !    call splev (knots_i, coefs_i, nknots_i, k, x, yhat_i, ext)
        ! end if
        !
        ! !$omp critical
        ! knots(:, j) = knots_i
        ! coefs(:, j) = coefs_i
        ! status(j) = status_i
        ! nknots(j) = nknots_i
        ! ssr(j) = ssr_i
        !  $omp end critical
    end do
    !$omp end do
    deallocate (ws)
    !$omp end parallel

    do i = 1, ns
       call print_report (iopt, svec(i), k, ssr(i), status(i), &
           nknots(i), knots(:, i), coefs(:, k), x, y, yhat(:, i), i)
    end do

end subroutine

subroutine print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, yhat, counter)
    integer, intent(in) :: iopt, k, n, counter
    integer (NF_ENUM_KIND), intent(in) :: status
    real (PREC), intent(in) :: s, ssr, knots(:), coefs(:)
    real (PREC), intent(in), dimension(:) :: x, y, yhat

    integer :: i, nstatus
    integer, dimension(bit_size(status)) :: istatus

    call decode_status (status, istatus, nstatus)

    print "(/,'(', i0, ')', t6, 'iopt: ', i2, '; smoothing factor: ', f6.1, '; spline degree: ', i1)", &
        counter, iopt, s, k
    print "(t6, 'status code: ', (i0, :, ','))", istatus(1:nstatus)
    ! Approximation is returned even if status /= NF_STATUS_OK as long as
    ! input arguments were valid.
    if (status /= NF_STATUS_INVALID_ARG) then
        print "(t6, 'SSR: ', es12.5e2)", ssr
        print "(t6, 'Number of knots: ', i0)", n
        print "(t6, 'Knots: ', *(t14, 10(f10.1, :, ', '), :, /))", knots(1:n)
        print "(t6, 'Coefs: ', *(t14, 10(g10.2e3, :, ', '), :, /))", coefs(1:n)
        print "(t6, a)", 'Evaluated points:'
        print "(t6, 5(3(a5, tr1), tr2))", ('x(i)', 'y(i)', 'sp(i)', i=1,5)
        print "(*(t6, 5(3(f5.1, :, tr1), tr2), :, /))", (x(i), y(i), yhat(i), i=1,size(x))
    end if
end subroutine

end program
