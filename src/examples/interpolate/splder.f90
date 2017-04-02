! Sample code demonstrating the use of splder to obtain derivates
! of fitted spline.
! Author: Richard Foltyn

program splder_demo

    use iso_fortran_env
    use numfort_interpolate, workspace => workspace_real64

    implicit none

    integer, parameter :: PREC = real64

    call example1 ()
contains

! Replicates the sample code in FITPACK's original mnspde.f file.
subroutine example1 ()

    integer, parameter :: m = 7, k_max = 5
    integer :: ext, k, i, n, nu
    type (status_t) :: status

    real (PREC), dimension(m) :: x
    real (PREC), dimension(m, k_max + 1) :: y
    real (PREC), dimension(:), allocatable :: knots, coefs
    type (workspace) :: ws
    real (PREC) :: ai

    ext = NF_INTERP_EVAL_EXTRAPOLATE

    ! initialize x
    x(1) = 0.0_PREC
    ai = 0.5d-1
    do i = 2, m
        x(i) = x(i-1) + ai
        ai = ai + 0.5d-1
    end do
    x(m) = 1_PREC

    ! loop through different spline degrees
    do k = 1, k_max
        ! number of knots
        n = 2 * (k + 1) + 4
        allocate (knots(n), coefs(n))
        ! initialize knots
        knots(1:k+1) = 0.0_PREC
        knots(k+6:) = 1_PREC
        ! interior knots
        knots(k+2) = 0.1_PREC
        knots(k+3) = 0.3_PREC
        knots(k+4) = 0.4_PREC
        knots(k+5) = 0.8_PREC

        ! generate b-spline coefficients
        ! This is weird as it leaves parts of coefs
        ! uninitialized, but that's how it's in the original
        ! FITPACK code
        do i = 1, n-k-1
            ai = real(i, PREC)
            coefs(i) = 0.1d-1 * ai * (ai - 5_PREC)
        end do

        ! evaluate spline derivates
        do i = 1, k+1
            nu = i - 1
            call splder (knots, coefs, n, k, nu, x, y(:, i), &
                ext=ext, work=ws, status=status)
        end do

        call print_report (n, k, knots, coefs, x, y, status)

        deallocate (knots, coefs)

    end do

end subroutine

subroutine print_report (n, k, knots, coefs, x, y, status)
    integer :: n, k
    real (PREC), dimension(:) :: knots, coefs, x
    real (PREC), dimension(:,:) :: y
    type (status_t), intent(in) :: status

    intent(in) :: k, knots, coefs, x, y

    integer, save :: ii = 1
    integer :: i, j

    integer, dimension(NF_MAX_STATUS_CODES) :: istatus
    integer :: nstatus

    call status_decode (status, istatus, nstatus)

    print "(/,'(', i0, ')', t6, 'spline degree: ', i1, '; error flag: ', i3)", &
        ii, k
    print "(t6, 'status code: ', *(i0, :, ', '))", istatus(1:nstatus)

    if (NF_STATUS_OK .in. status) then
        print "(t6, 'Number of knots: ', i0)", n
        print "(t6, 'Knots: ', *(t14, 10(f8.1, :, ', '), :, /))", knots
        print "(t6, 'Coefs: ', *(t14, 10(f8.4, :, ', '), :, /))", coefs
        print "(t6, a)", 'Evaluated derivatives:'
        print "(t8, a4, tr2, *(tr4, 'sp_', i1, '(i)', :, tr2))", 'x(i)', (i, i=0,k)

        do i = 1, size(x)
            print "(t8, f4.2, tr2, *(e11.4e2, :, tr2))", x(i), (y(i, j), j=1,k+1)
        end do
    end if

    ii = ii + 1
end subroutine

end program
