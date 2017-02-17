! Demo code for wrapper for FITPACK's concon routine to fit cubic splines
! with convexity/concavity constraints

program demo_concon

    use iso_fortran_env
    use numfort_interpolate

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()

contains

! EXAMPLE 1 replicates the code contained in mncoco.f in FITPACK.
! This fits a function that is concave in all intervals for various smoothing
! parameters.
subroutine example1 ()

    ! number of points, degree of spline (cubic)
    integer, parameter :: m = 16, k = 3

    real (PREC), dimension(m) :: x, y, w, v, sx, sx2, s1, s2
    real (PREC), dimension(:), allocatable :: knots, coefs
    real (PREC) :: ssr
    logical, dimension(:), allocatable :: bind
    integer :: status, iopt, maxtr, maxbin, is, nest, n, status_eval
    ! different smoothness parameters
    real (PREC) :: s(3) = [real(PREC) :: 0.2, 0.04, 0.0002]

    x = [real(PREC) :: 0.1, 0.3, 0.5, 0.7, 0.9, 1.25, 1.75, 2.25, 2.75, 3.5, &
        4.5, 5.5, 6.5, 7.5, 8.5, 9.5]

    y = [real(PREC) :: 0.124,0.234,0.256,0.277,0.278, 0.291, 0.308, 0.311, &
        0.315, 0.322, 0.317, 0.326, 0.323, 0.321, 0.322, 0.328]

    w = 1.0
    w(1) = 10
    w(2) = 3
    w(16) = 10

    ! "globally" concave approximation
    v = 1.0

    iopt = 0
    maxtr = 100
    maxbin = 10

    nest = concon_get_nest (m)
    allocate (knots(nest), coefs(nest), bind(nest))

    do is = 1, size(s)
        call concon (x, y, v, s(is), n, knots, coefs, iopt, w, maxtr, maxbin, &
            sx=sx, bind=bind, ssr=ssr, status=status)

        ! evaluate spline using splev to check values in sx
        call splev (knots, coefs, n, 3, x, sx2, status=status)
        if (status /= INTERP_STATUS_SUCCESS) then
            write (ERROR_UNIT, *) "Failed to evaluate spline"
        end if
        if (any(abs(sx-sx2) > 1d-10)) then
            write (ERROR_UNIT, *) "s(x) returned from CONCON and SPLEV differ"
        end if

        ! evaluate first- and second-order derivatives
        call splder (knots(1:n), coefs(1:n), k=k, order=1, x=x, y=s1, &
            ext=INTERP_EVAL_ERROR, status=status_eval)

        call splder (knots, coefs, n, k, 2, x, s2, &
            INTERP_EVAL_ERROR, status=status_eval)

        call print_report (iopt, s(is), ssr, status, n, knots, coefs, x, y, &
            sx, s1, s2, bind)
    end do

end subroutine

subroutine print_report (iopt, s, ssr, status, n, knots, coefs, x, y, yhat, &
        s1, s2, bind, counter)
    integer, intent(in) :: iopt, status, n, counter
    real (PREC), intent(in) :: s, ssr, knots(:), coefs(:)
    real (PREC), intent(in), dimension(:) :: x, y, yhat, s1, s2
    logical, dimension(:), intent(in) :: bind

    optional :: counter
    integer, save :: ii = 1
    integer :: i
    character (len=size(bind)) :: str_bind

    if (present(counter)) ii = counter

    ! determine knots where convexity restriction is binding
    str_bind = ''
    do i = 1, size(knots)-6
        if (bind(i)) str_bind(i+3:i+3) = '*'
    end do

    print "(/,'(', i0, ')', t6, 'iopt: ', i2, '; smoothing factor: ', es10.1e2)", &
        ii, iopt, s
    print "(t6, 'SSR: ', es15.5e3, '; error flag: ', i3)", ssr, status
    print "(t6, 'Number of knots: ', i0)", n
    ! print knots where convexitiy restriction is binding with an adjacent *
    print "(t6, 'Knots: ', *(t14, 10(f8.1, :, a1, tr1), :, /))", &
        (knots(i), str_bind(i:i), i=1,size(knots))
    print "(t6, 'Coefs: ', *(t14, 10(f8.4, :, tr2), :, /))", coefs(1:n)
    print "(t6, a)", 'Evaluated points:'
    print "(t6, 2(3(a5, tr1), 2(a10, tr1), :, tr2))", &
        ('x(i)', 'y(i)', 's(i)', 's1(i)', 's2(i)', i=1,2)
    print "(*(t6, 2(3(f5.1, tr1), 2(es10.2e3, tr1), :, tr2), :, /))", &
        (x(i), y(i), yhat(i), s1(i), s2(i), i=1,size(x))

    ii = ii + 1
end subroutine

end program
