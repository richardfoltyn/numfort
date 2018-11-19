


subroutine __APPEND(hermegauss_check_input,__PREC) (x, w, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: w
    type (status_t), intent(out) :: status

    if (size(x) == 0) goto 100
    if (size(w) /= size(x)) goto 100

    ! No invalid input encountered
    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG
end subroutine



pure subroutine __APPEND(hermecompanion,__PREC) (c, mat)
    !*  HERMECOMPANION returns the scaled companion matrix of C.
    !
    !   The basis polynomials are scaled do that the companion matrix is
    !   symmetric when C is an Hermite_e basis polynomial. The provides
    !   better eigenvalue estimates than the unscaled case and for basis
    !   polynomials the eigenvalues are guaranteed to be real if
    !   LAPACK's SYEVD is used to obtain them.
    !
    !   Note: this is a Fortran port of the Python implementation in
    !   Numpy's hermecompanion() function.
    !
    !   Note: this is an internal routine, so no input validation is performed.
    !
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: c
        !*  Array of Hermite_e series coefficients ordered from low to high
        !   degree.
    real (PREC), intent(out), dimension(:,:) :: mat
        !*  Scaled companion matrix of dimension (deg, deg) where
        !   deg = SIZE(c) - 1.

    integer :: n, i
    real (PREC), dimension(:), allocatable :: rwork1d
    real (PREC) :: x

    n = size(c) - 1

    ! Abuse first column of MAT as workspace
    ! Perform two steps in Numpy code in one go, ie creating the values
    ! 1/sqrt(i) and computing the cumulative product (but not yet reversing
    ! the element order)
    mat(1,1) = 1.0
    do i = 2, n
        mat(i,1) = mat(i-1,1) / sqrt(real(n-i+1, PREC))
    end do

    allocate (rwork1d(n))
    do i = 1, n
        rwork1d(i) = mat(n-i+1,1) * c(i) / c(n+1)
    end do

    mat = 0.0
    do i = 1, n-1
        x = sqrt(real(i, PREC))
        ! Upper "diagonal"
        mat(i,i+1) = x
        ! Lower "diagonal"
        mat(i+1,i) = x
    end do

    ! Adjust last column
    do i = 1, n
        mat(i,n) = mat(i,n) - rwork1d(i)
    end do

    deallocate (rwork1d)

end subroutine




pure subroutine __APPEND(eval_normed_herme,__PREC) (x, deg, hx)
    !*  EVAL_NORMED_HERME evaluates a normalized Hermite_e polynomial
    !   of degree N at a set of given points X.
    !
    !   Note: this is a port of the Numpy function _normed_hermite_e_n
    !   to Fortran.
    !
    !   This function is meant for internal use only, hence no input validation
    !   is performed.
    integer, parameter :: PREC = __PREC
    real (PREC), dimension(:), intent(in), contiguous :: x
        !*  Points at which to evaluate polynomial
    integer, intent(in) :: deg
        !*  Degree of the normalized Hermite_e polynomial to be evaluated.
    real (PREC), dimension(:), intent(out), contiguous :: hx
        !*  Polynomial values, assumed to be the same shape as X

    real (PREC), parameter :: PI = __APPEND(PI,__PREC)
    real (PREC), dimension(:), allocatable :: c0, c1
    real (PREC) :: v
    real (PREC) :: rdeg
    integer :: i, n

    ! Numpy's 1/sqrt(sqrt(2*pi))
    v = (2.0_PREC * PI)**(-1.0_PREC/4.0_PREC)

    if (deg == 0) then
        hx = v
        return
    end if

    n = size(x)
    allocate (c0(n), source=0.0_PREC)
    allocate (c1(n), source=v)

    rdeg = real(deg, PREC)

    do i = 1, deg - 1
        ! Abuse HX as temporary array
        hx = c0
        c0(:) =  - c1 * sqrt((rdeg - 1.0_PREC)/rdeg)
        c1(:) = hx + c1 * x * sqrt(1.0_PREC/rdeg)
        rdeg = rdeg - 1.0
    end do

    hx = c0 + c1 * x

end subroutine



subroutine __APPEND(hermegauss,__PREC) (x, w, status)
    !*  HERMEGAUSS computes the sample points and weights for Gauss-Hermite_e
    !   quadrature. These sample points and weights will correctly integrate
    !   polynomials of degree 2*n-1 or less, where n=SIZE(X), over the
    !   interval (-inf, inf) with the weight function exp(-x^2/2).
    !
    !   Note that the weight function is the standard normal density
    !   scaled by sqrt(2*pi), thus the weights have to be divided by
    !   this value to compute expectations of (functions of) standard
    !   normal random variables.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(out), dimension(:), contiguous :: x
        !*  Array containing the sample points
    real (PREC), intent(out), dimension(:), contiguous :: w
        !*  Array containing the weights
    type (status_t), intent(out), optional :: status
        !*  Optional exit status

    type (status_t) :: lstatus
    real (PREC), dimension(:,:), allocatable :: mat
    real (PREC), dimension(:), allocatable :: rwork1d, dy, df
    integer :: i, n, deg
    real (PREC) :: wmax, val
    real (PREC), parameter :: PI = __APPEND(PI,__PREC)
    ! Arguments used only for all to SYEVD
    integer, dimension(:), allocatable :: iwork1d
    integer :: liwork, lwork, info

    call hermegauss_check_input (x, w, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    deg = size(x)

    allocate (rwork1d(deg+1), source=0.0_PREC)
    rwork1d(deg+1) = 1.0_PREC

    ! Compute scaled companion matrix
    allocate (mat(deg,deg))
    call hermecompanion (rwork1d, mat)
    deallocate (rwork1d)

    ! Compute eigenvalues of symmetric matrix MAT
    ! Step one: workspace query
    allocate (rwork1d(1))
    allocate (iwork1d(1))
    call LAPACK_SYEVD (jobz='N', uplo='L', n=deg, a=mat, lda=deg, w=x, &
        work=rwork1d, lwork=-1, iwork=iwork1d, liwork=-1, info=info)

    lwork = int(rwork1d(1))
    liwork = iwork1d(1)
    deallocate (rwork1d, iwork1d)
    allocate (rwork1d(lwork))
    allocate (iwork1d(liwork))
    ! Actuall call to compute (only!) eigenvalues of MAT
    call LAPACK_SYEVD (jobz='N', uplo='L', n=deg, a=mat, lda=deg, w=x, &
        work=rwork1d, lwork=lwork, iwork=iwork1d, liwork=liwork, info=info)

    deallocate (rwork1d, iwork1d)

    ! Improve roots by one application of Newton
    allocate (df(deg), dy(deg))
    call eval_normed_herme (x, deg, dy)
    call eval_normed_herme (x, deg - 1, df)
    x = x - dy / df / sqrt(real(deg, PREC))
    deallocate (df, dy)

    ! compute the weights: We scale the factor to avoid possible numerical
    ! overflow.
    call eval_normed_herme (x, deg - 1, w)
    wmax = maxval(abs(w))
    w = w / wmax
    w = 1.0_PREC / w**2.0_PREC

    ! Symmetrize weights and x's
    n = int((deg+1)/2)
    do i = 1, n
        val = (w(i) + w(deg-i+1)) / 2.0_PREC
        w(i) = val
        w(deg-i+1) = val

        val = (x(i) - x(deg-i+1)) / 2.0_PREC
        x(deg-i+1) = - val
        x(i) = val
    end do

    ! Scale up weights since we are computing quadrature data for the
    ! non-normalized weighting function exp(-x**2/2)
    w = w * (sqrt(2.0_PREC * PI) / sum(w))

100 continue
    if (present(status)) status = lstatus

end subroutine
