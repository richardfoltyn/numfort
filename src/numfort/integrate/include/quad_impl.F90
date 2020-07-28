

subroutine quad_check_input (a, b, epsabs, epsrel, nmax, status)
    !*  Input checker routine for QUAD.

    real (PREC), intent(in) :: a, b
        !*  Integration bounds. Note that QUADPACK does not require a <= b.
    real (PREC), intent(in) :: epsabs, epsrel
    integer, intent(in) :: nmax
    type (status_t), intent(out) :: status

    status = NF_STATUS_INVALID_ARG

    if (epsabs < 0.0_PREC .or. epsrel < 0.0_PREC) goto 100
    if (nmax < 1) goto 100

    status = NF_STATUS_OK

100 continue

end subroutine


subroutine quad (fcn, a, b, res, epsabs, epsrel, nmax, work, &
        abserr, nfev, status)
    !*  Integrate function on interval [a,b] using QUADPACK routines.

    procedure (f_integrand) :: fcn
        !*  Function defining the integrand
    real (PREC), intent(in) :: a
        !*  Lower bound of integration
    real (PREC), intent(in) :: b
        !*  Upper bound of integration
    real (PREC), intent(out) :: res
        !*  Approximated integral
    real (PREC), intent(in), optional :: epsabs
        !*  Absolute error tolerance
    real (PREC), intent(in), optional :: epsrel
        !*  Relative error tolerance
    integer, intent(in), optional :: nmax
        !*  Upper bound on the number of subintervals used in adaptive algorithm
    type (workspace), intent(inout), optional :: work
        !*  Workspace object (optional). If not present, routine will automatically
        !   allocate required arrays.
    real (PREC), intent(out), optional :: abserr
        !*  If present, contains estimate of the absolute error in the result on exit.
    integer, intent(out), optional :: nfev
        !*  If present, contains number of integrand evaluations on exit.
    type (status_t), intent(out), optional :: status
        !*  If present, contains exit status.

    type (workspace), pointer :: ptr_work
    real (PREC) :: lepsabs, lepsrel
    integer :: lnmax

    type (status_t) :: lstatus
    integer :: niwrk, nrwrk
    ! Arguments to QUADPACK routines
    integer :: ier, neval
    integer :: ln
        !*  Actual number of subintervals used
    real (PREC) :: labserr
    real (PREC), dimension(:), pointer, contiguous :: ptr_alist, ptr_blist, &
        ptr_rlist, ptr_elist
    integer, dimension(:), pointer, contiguous :: ptr_iord

    nullify (ptr_work)
    nullify (ptr_alist, ptr_blist, ptr_rlist, ptr_elist, ptr_iord)

    res = 0.0_PREC
    labserr = 0.0_PREC
    lstatus = NF_STATUS_OK

    ! Use Scipy's tolerance levels by default
    lepsabs = sqrt(epsilon(0.0_PREC))
    lepsrel = sqrt(epsilon(0.0_PREC))
    lnmax = 50

    if (present(epsabs)) lepsabs = epsabs
    if (present(epsrel)) lepsrel = epsrel
    if (present(nmax)) lnmax = nmax

    call quad_check_input (a, b, lepsabs, lepsrel, lnmax, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    nrwrk = lnmax * 4
    niwrk = lnmax

    ! Allocate workspace arrays
    call assert_alloc_ptr (work, ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk, niwrk=niwrk)

    ptr_alist => ptr_work%rwrk(1:lnmax)
    ptr_blist => ptr_work%rwrk(lnmax+1:2*lnmax)
    ptr_rlist => ptr_work%rwrk(2*lnmax+1:3*lnmax)
    ptr_elist => ptr_work%rwrk(3*lnmax+1:4*lnmax)
    ptr_iord => ptr_work%iwrk(1:lnmax)

    ! Call underlying QUADPACK implementation
    call qagse (fcn, a, b, lepsabs, lepsrel, lnmax, res, labserr, neval, ier, &
        ptr_alist, ptr_blist, ptr_rlist, ptr_elist, ptr_iord, ln)

    call map_qagse_status (ier, lstatus)


100 continue

    ! Clean up any "local" working arrays
    call assert_dealloc_ptr (work, ptr_work)

    ! Set optional return values
    if (present(status)) status = lstatus
    if (present(abserr)) abserr = labserr
    if (present(nfev)) nfev = neval

end subroutine



pure subroutine map_qagse_status (ier, status)
    !*  MAP_QAGSE_STATUS maps the integer value returned by
    !   QUADPACK's QAGSE routine to NUMFORT status codes.

    integer, intent(in) :: ier
    type (status_t), intent(inout) :: status

    select case (ier)
    case (0)
        status = NF_STATUS_OK
    case (1)
        status = NF_STATUS_MAX_INTERVALS
    case (2)
        status = NF_STATUS_NOT_CONVERGED
    case (3)
        status = NF_STATUS_INTEGRAND_ERROR
    case (4)
        status = NF_STATUS_NOT_CONVERGED
    case (5)
        status = NF_STATUS_INTEGRAND_ERROR
    case (6)
        status = NF_STATUS_INVALID_ARG
    end select

    status%code_orig = ier
end subroutine
