! Sample code demonstrating the use of concon, splev and splder from
! the numfort_interpolate module, using OpenMP.
! Author: Richard Foltyn

program curfit_demo_omp

    use iso_fortran_env
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

    real (PREC), dimension(m) :: x, y, eps, w, v
    real (PREC), dimension(:,:), allocatable :: sp, sp1, sp2
    real (PREC), dimension(ns) :: svec
    integer :: i, nest, iopt, nseed, ext
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
    real (PREC), dimension(:,:), allocatable :: sp_2, sp1_2, sp2_2

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
    iopt = 0
    w = 1/abs(eps)
    ! evaluate spline at original x points
    ext = NF_INTERP_EVAL_BOUNDARY
    ! fit concave spline
    v = 1.0_PREC

    nest = concon_get_nest (m)
    allocate (knots(nest, ns), coefs(nest, ns))
    allocate (sp(m, ns), sp1(m, ns), sp2(m, ns))

    !$omp parallel default(none) private(j, ws, s, tid) &
    !$omp shared(knots2, coefs2, nknots2, ssr2, status2) &
    !$omp shared(svec, x, y, v, knots, coefs, nknots, iopt, w, ssr, status, nest) &
    !$omp shared(sp, sp1, sp2, ext)
    tid = omp_get_thread_num ()
    allocate (ws)

    !$omp do
    do j = 1, ns
        s = svec(j)
        call fit (x, y, v, s, knots(:, j), coefs(:, j), nknots(j), &
            ws, ssr(j), status(j), sp(:, j), sp1(:, j), sp2(:, j), &
            iopt, w, ext=ext)
    end do
    !$omp end do
    deallocate (ws)
    !$omp end parallel

    do i = 1, ns
       call print_report (iopt, svec(i), ssr(i), status(i), &
           nknots(i), knots(:, i), coefs(:, i), x, y, sp(:, i), sp1(:, i), &
           sp2(:, i), i)
    end do

    ! compute sequentially
    allocate (knots2(nest, ns), coefs2(nest, ns))
    allocate (sp_2(m, ns), sp1_2(m, ns), sp2_2(m, ns))
    allocate (ws)

    do j = 1, ns
        s = svec(j)
        call fit (x, y, v, s, knots2(:, j), coefs2(:, j), nknots2(j), &
            ws, ssr2(j), status2(j), sp_2(:, j), sp1_2(:, j), sp2_2(:, j), &
            iopt, w, ext=ext)
    end do

    deallocate (ws)

    print *, "== Summary of OpenMP vs. sequential results =="
    do j = 1, ns
        print '(i3, a, es9.2e2, "; ", a, es9.2e2, "; ", a, l1, "; ", a, es9.2e2, "; ", a, l1)', &
            j, &
            "d(knots)=", maxval(abs(knots(:,j)-knots2(:,j))), &
            "d(coefs)=", maxval(abs(coefs(:,j)-coefs2(:,j))), &
            "knots. eq=", nknots(j) == nknots2(j), &
            "d(ssr)=", abs(ssr(j)-ssr2(j)), &
            "status eq.=", status(j) == status2(j)
    end do

end subroutine

pure subroutine fit (x, y, v, s, knots, coefs, nknots, ws, ssr, status, &
        sp, sp1, sp2, iopt, w, ext)
    real (PREC), intent(in), dimension(:) :: x, y, v
    real (PREC), intent(in) :: s
    real (PREC), intent(out), dimension(:) :: knots, coefs
    integer, intent(out) :: nknots
    class (workspace), intent(in out) :: ws
    real (PREC), intent(out) :: ssr
    type (status_t), intent(out) :: status
    real (PREC), intent(out), dimension(:) :: sp, sp1, sp2
    integer, intent(in), optional :: iopt
    real (PREC), intent(in), dimension(:), optional :: w
    integer (NF_ENUM_KIND), intent(in), optional :: ext

    logical :: stat_ok

    call concon (x, y, v, s, knots, coefs, nknots, &
       iopt=iopt, work=ws, w=w, ssr=ssr, sx=sp, status=status)

    stat_ok = any ([NF_STATUS_OK, NF_STATUS_APPROX] .in. status)

    if (stat_ok) then
       call splder (knots, coefs, nknots, 3, 1, x, sp1, work=ws, ext=ext)
       call splder (knots, coefs, nknots, 3, 2, x, sp2, work=ws, ext=ext)
    end if
end subroutine

subroutine print_report (iopt, s, ssr, status, n, knots, coefs, x, y, &
        sp, sp1, sp2, counter)
    integer, intent(in) :: iopt, n, counter
    type (status_t), intent(in) :: status
    real (PREC), intent(in) :: s, ssr, knots(:), coefs(:)
    real (PREC), intent(in), dimension(:) :: x, y, sp, sp1, sp2
    logical :: stat_ok

    integer :: i, nstatus
    integer, dimension(NF_MAX_STATUS_CODES) :: istatus

    call status%decode (istatus, nstatus)

    stat_ok = any ([NF_STATUS_OK, NF_STATUS_APPROX] .in. status)
    
    print "(/,'(', i0, ')', t6, 'iopt: ', i2, '; smoothing factor: ', es12.5e2)", &
        counter, iopt, s
    print "(t6, 'status code: ', *(i0, :, ', '))", istatus(1:nstatus)
    ! Approximation is returned even if status /= NF_STATUS_OK as long as
    ! input arguments were valid.
    if (stat_ok) then
        print "(t6, 'SSR: ', es12.5e2)", ssr
        print "(t6, 'Number of knots: ', i0)", n
        print "(t6, 'Knots: ', *(t14, 8(f8.3, :, ', '), :, /))", knots(1:n)
        print "(t6, 'Coefs: ', *(t14, 8(f8.3, :, ', '), :, /))", coefs(1:n)
        print "(t6, a)", 'Evaluated points:'
        print "(t6, 2(5(a7, tr1), tr5))", ('x(i)', 'y(i)', 'sp(i)', 'sp1(i)', 'sp2(i)', i=1,2)
        print "(*(t6, 2(5(f7.2, :, tr1), tr5), :, /))", &
            (x(i), y(i), sp(i), sp1(i), sp2(i), i=1,size(x))
    else
        print '(t6, a)', "No spline approximation returned"
    end if
end subroutine

end program
