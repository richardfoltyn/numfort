! Sample code demonstrating the use of curfit and splev from
! the numfort_interpolate module.
! Author: Richard Foltyn

program curfit_demo

    use iso_fortran_env
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
    integer :: i, k, nest_max, n, status, j
    real (PREC) :: ssr
    real (PREC), dimension(:), allocatable :: knots, coefs
    type (workspace) :: ws

    ! initialize data
    x = [ real(PREC) :: (i, i = 0, m-1)]
    y = [ real(PREC) :: 1.0,1.0,1.4,1.1,1.0,1.0,4.0,9.0,13.0,13.4,12.8,13.1,13.0, &
        14.0,13.0,13.5,10.0,2.0,3.0,2.5,2.5,2.5,3.0,4.0,3.5]

    ! max. number of elements needed for knots/coefs
    nest_max = curfit_get_nest(m=m, k=5)
    allocate (knots(nest_max), coefs(nest_max))

    print "(/, '### Example ', i0)", 1

    do k = 3,5,2
        ! iterate through various smoothing parameters and interpolation tasks
        do i = 1, size(svec)
            s = svec(i)
            iopt = ioptvec(i)

            call curfit (x, y, k, s, n, knots, coefs, iopt, ssr=ssr, work=ws, &
                status=status)

            ! test different call syntax without specifying n directly
            call splev (knots(1:n), coefs(1:n), k=k, x=x, y=yhat, &
                ext=INTERP_EVAL_EXTRAPOLATE, status=status)
            call print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, yhat)
        end do

        ! least-squares fitting
        s = 0.0_PREC
        iopt = -1
        knots(k+2:k+8) = [(3.0_PREC * j, j = 1,7)]
        n = 9 + 2*k
        call curfit (x, y, k, s, n, knots, coefs, iopt, ssr=ssr, &
            work=ws, status=status)

        call splev (knots, coefs, n, k, x, yhat, INTERP_EVAL_EXTRAPOLATE, status)
        call print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, yhat)
    end do

end subroutine

subroutine example2 ()

    ! number of points
    integer, parameter :: m = 25

    real (PREC), dimension(m) :: x, y, yhat, eps, w
    real (PREC) :: s
    integer :: i, k, nest, n, status, iopt, nseed, ext
    real (PREC) :: ssr
    real (PREC), dimension(:), allocatable :: knots, coefs, seed
    type (workspace) :: ws

    print "(/, '### Example ', i0)", 2

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
    s = 100.0_PREC
    w = 1/eps

    nest = curfit_get_nest (m=m, k=k)
    allocate (knots(nest), coefs(nest))

    call curfit (x, y, k, s, n, knots, coefs, iopt, w=w, work=ws, ssr=ssr, &
        status=status)

    ! evaluate spline at original x points
    ext = INTERP_EVAL_BOUNDARY
    call splev (knots, coefs, n, k, x, yhat, ext, status)

    call print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, yhat, &
        counter=1)

end subroutine

subroutine print_report (iopt, s, k, ssr, status, n, knots, coefs, x, y, yhat, &
        counter)
    integer, intent(in) :: iopt, k, status, n, counter
    real (PREC), intent(in) :: s, ssr, knots(:), coefs(:)
    real (PREC), intent(in), dimension(:) :: x, y, yhat

    optional :: counter
    integer, save :: ii = 1
    integer :: i

    if (present(counter)) ii = counter

    print "(/,'(', i0, ')', t6, 'iopt: ', i2, '; smoothing factor: ', f6.1, '; spline degree: ', i1)", &
        ii, iopt, s, k
    print "(t6, 'SSR: ', es12.5e2, '; error flag: ', i3)", ssr, status
    print "(t6, 'Number of knots: ', i0)", n
    print "(t6, 'Knots: ', *(t14, 10(f8.1, :, ', '), :, /))", knots(1:n)
    print "(t6, 'Coefs: ', *(t14, 10(f8.4, :, ', '), :, /))", coefs(1:n)
    print "(t6, a)", 'Evaluated points:'
    print "(t6, 5(3(a5, tr1), tr2))", ('x(i)', 'y(i)', 'sp(i)', i=1,5)
    print "(*(t6, 5(3(f5.1, :, tr1), tr2), :, /))", (x(i), y(i), yhat(i), i=1,size(x))

    ii = ii + 1
end subroutine

end program
