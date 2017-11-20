real (PREC), intent(out), dimension(:) :: x
    !*  Array to store evenly spaced numbers over specified interval.
real (PREC), intent(in) :: xfrom
    !*  Starting value of sequence.
real (PREC), intent(in) :: xto
    !*  End value of sequence to be created. If step size is explicitly given
    !   and the array x is not large enough, the last element of X can be
    !   significantly different from XTO.
real (PREC), intent(in), optional :: step
    !*  Optional step size used to create evenly spaced elements. If not present,
    !   step size is imputed from the boundary values and the size of argument X.
integer, intent(out), optional :: res_n
    !*  If present, contains the actual number of elements returned. If
    !   the STEP argument is not present, res_n = size(x).
real (PREC), intent(out), optional :: res_step
    !*  If present, contains the actual step size used. If the STEP argument
    !   is present, then res_step = step.

real (PREC) :: lstep, xnext, sgn
integer :: i, nx, n

n = 0
nx = size(x)
lstep = 0.0_PREC

if (nx == 0) goto 100

if (.not. present(step)) then
    ! Infer step size from number boundaries and array size.
    lstep = (xto - xfrom) / (nx - 1.0_PREC)

    x(1) = xfrom
    do i = 2, nx-1
        x(i) =  xfrom + (i-1) * lstep
    end do

    ! avoid rounding errors at end point
    x(nx) = xto
    ! Actual number of elements used is identical to array size
    n = nx
else
    sgn = signum (xto-xfrom)
    lstep = step

    xnext = xfrom
    do i = 1, nx
        x(i) = xnext

        xnext = xnext + lstep
        if (xnext*sgn > xto*sgn) exit
    end do

    ! Actual number of elements used need not be identical to array size
    n = min(i, nx)
end if

100 continue
    if (present(res_n)) res_n = n
    if (present(res_step)) res_step = lstep
