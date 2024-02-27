! Sample code demonstrating the use of curfit and splev from
! the numfort_interpolate module.
! Author: Richard Foltyn

program curfit_demo

    use iso_fortran_env
    use numfort_common, workspace => workspace_real64
    use numfort_interpolate

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()
    call example2 ()

contains

! replicates original fitpack examples from mncurf.f, see this file
! for some hints of what is going on exactly.
! For a given set of (x,y) points, fit splines of degree 3 and 5 using
! various smoothness parameters, as well as LS.
! Report results and fitted points.
subroutine example1 ()

    ! number of points
    integer, parameter :: m = 25

    real (PREC), dimension(m) :: x, y, yhat
    real (PREC) :: s, svec(6) = [ real(PREC) :: 1000, 60, 10, 30, 30, 0]
    integer :: iopt, ioptvec(6) = [0, 1, 1, 1, 0, 0]
    integer :: i, k, nest_max, n, j, ext
    type (status_t) :: status
    real (PREC) :: ssr
    real (PREC), dimension(:), allocatable :: knots, coefs
    type (workspace) :: ws

    ! initialize data
    x = [ real(PREC) :: (i, i = 0, m-1)]
    y = [ real(PREC) :: 1.0,1.0,1.4,1.1,1.0,1.0,4.0,9.0,13.0,13.4,12.8,13.1,13.0, &
        14.0,13.0,13.5,10.0,2.0,3.0,2.5,2.5,2.5,3.0,4.0,3.5]

    ext = NF_INTERP_EVAL_EXTRAPOLATE
    ! max. number of elements needed for knots/coefs
    nest_max = curfit_get_nest(m=m, k=5)
    allocate (knots(nest_max), coefs(nest_max))

    print "(/, '### Example ', i0)", 1

    do k = 3,5,2
        ! iterate through various smoothing parameters and interpolation tasks
        do i = 1, size(svec)
            s = svec(i)
            iopt = ioptvec(i)

            call curfit (x, y, k, s, knots, coefs, n, iopt=iopt, ssr=ssr, &
                work=ws, status=status)

            ! test different call syntax without specifying n directly
            call splev (knots(1:n), coefs(1:n), k=k, x=x, y=yhat, &
                ext=ext, status=status)
            call print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, yhat)
        end do

        ! least-squares fitting
        s = 0.0_PREC
        iopt = -1
        knots(k+2:k+8) = [(3.0_PREC * j, j = 1,7)]
        n = 9 + 2*k
        call curfit (x, y, k, s, knots, coefs, n, iopt, ssr=ssr, &
            work=ws, status=status)

        call splev (knots(1:n), coefs(1:n), k, x, yhat, ext, status)
        call print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, yhat)
    end do

end subroutine

subroutine example2 ()

    ! number of points
    integer, parameter :: m = 25
    integer, parameter :: ns = 10

    real (PREC), dimension(m) :: x, y, yhat, eps, w
    real (PREC), dimension(10) :: svec
    integer :: i, k, nest, n, iopt, nseed, ext
    type (status_t) :: status
    real (PREC) :: ssr, s
    real (PREC), dimension(:), allocatable :: knots, coefs, seed
    type (workspace), pointer :: ws
    logical :: stat_ok

    print "(/, '### Example ', i0)", 2

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
    w = 1/abs(eps + 1.0d-5)
    ! evaluate spline at original x points
    ext = NF_INTERP_EVAL_BOUNDARY

    nest = curfit_get_nest (m=m, k=k)
    allocate (knots(nest), coefs(nest))

    do i = 1, ns
        s = svec(i)
        allocate (ws)

        call curfit (x, y, k, s, knots, coefs, n, w=w, iopt=iopt, &
            work=ws, ssr=ssr, status=status)

        stat_ok = .not. (NF_STATUS_INVALID_ARG .in. status)

        if (stat_ok) then
            call splev (knots(1:n), coefs(1:n), k, x, yhat, ext)
        end if

        if (i == 1) then
            call print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, &
                yhat, counter=1)
        else
            call print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, &
                yhat)
        end if

        deallocate (ws)
    end do


end subroutine

subroutine print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, yhat, &
        counter)
    integer, intent(in) :: iopt, k, n, counter
    type (status_t), intent(in) :: status
    real (PREC), intent(in) :: s, ssr, knots(:), coefs(:)
    real (PREC), intent(in), dimension(:) :: x, y, yhat

    optional :: counter
    integer, save :: ii = 1
    integer :: i, nstatus
    integer, dimension(NF_MAX_STATUS_CODES) :: istatus
    logical :: stat_ok

    if (present(counter)) ii = counter

    call status_decode (status, istatus, nstatus)
    stat_ok = .not. (NF_STATUS_INVALID_ARG .in. status)

    print "(/,'(', i0, ')', t6, 'iopt: ', i2, '; smoothing factor: ', es12.5e2, '; spline degree: ', i1)", &
        ii, iopt, s, k
    print "(t6, 'status codes: ', *(i0, :, ', '))", istatus(1:nstatus)
    if (stat_ok ) then
        print "(t6, 'SSR: ', es12.5e2)", ssr
        print "(t6, 'Number of knots: ', i0)", n
        print "(t6, 'Knots: ', *(t14, 8(f8.3, :, ', '), :, /))", knots(1:n)
        print "(t6, 'Coefs: ', *(t14, 8(f8.3, :, ', '), :, /))", coefs(1:n)
        print "(t6, a)", 'Evaluated points:'
        print "(t6, 5(3(a5, tr1), tr2))", ('x(i)', 'y(i)', 'sp(i)', i=1,5)
        print "(*(t6, 5(3(f5.1, :, tr1), tr2), :, /))", (x(i), y(i), yhat(i), i=1,size(x))
    end if

    ii = ii + 1
end subroutine

end program
