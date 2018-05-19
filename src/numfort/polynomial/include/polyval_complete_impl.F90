

pure subroutine __APPEND(polybasis_check_input,__PREC) (x, k, basis, status)
    !*  POLYBASIS_CHECK_INPUT performs input checks for user-provided
    !   arguments for POLYBASIS_COMPLETE.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: x
    integer, intent(in) :: k
    real (PREC), intent(out), dimension(:,:) :: basis
    type (status_t), intent(out), optional :: status

    integer :: npoints, nd, nterms

    status = NF_STATUS_OK

    npoints = size(x, 2)
    nd = size(x, 1)

    call check_nonneg (k, k, "k", status=status)
    if (status /= NF_STATUS_OK) goto 100

    nterms = poly_complete_get_nterms (nd, k)

    if (size(basis, 1) < nterms .or. size(basis, 2) < npoints) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

100 continue

end subroutine


pure subroutine __APPEND(polybasis,__PREC) (x, k, basis, status)
    !*  POLYBASIS_COMPLETE returns the basis functions (or terms) of a
    !   complete polynomial of given degree and number of variables,
    !   evaluated at a given set of x.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x
        !*  Set of x where basis functions should be evaluated
        !   (each column contains on point of dimension N)
    integer, intent(in) :: k
        !*  (Maximum) degree of complete polynomial
    real (PREC), intent(out), dimension(:,:) :: basis
        !*  Basis functions evaluated at X
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    integer, dimension(:,:), allocatable :: exponents
    real (PREC), dimension(:), allocatable :: xi
    integer :: nbasis, nd, i, j, nx, ib
    real (PREC) :: bi

    lstatus = NF_STATUS_OK

    call polybasis_check_input (x, k, basis, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! number of points and dimensions
    nd = size(x, 1)
    nx = size(x, 2)
    ! Number of basis functions
    nbasis = poly_complete_get_nterms (nd, k)

    ! No points given, hence nothing to do
    if (nx == 0) goto 100

    ! Basis functions for 0-dimensional data: in this case we should just
    ! return the constant term and exit immediately.
    if (nd == 0) then
        basis(1,:) = 1.0
        goto 100
    end if

    allocate (xi(nd))

    ! Array that contains all permutations of exponents that occur in
    ! given complete polynomial (one per column).
    allocate (exponents(nd, nbasis))

    ! obtain sorted permutations of exponents that sum to k=2,...,k
    call polyexponents_complete (nd, k, exponents, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    do j = 1, nx
        xi(:) = x(:, j)

        do ib = 1, nbasis
            bi = 1.0
            do i = 1, nd
                bi = bi * xi(i) ** exponents(i, ib)
            end do
            basis(ib, j) = bi
        end do
    end do


100 continue

    if (allocated(exponents)) deallocate (exponents)
    if (allocated(xi)) deallocate (xi)

    if (present(status)) status = lstatus

end subroutine


pure subroutine __APPEND(polybasis_scalar,__PREC) (x, k, basis, status)
    !*  POLYBASIS_COMPLETE_SCALAR returns the basis functions (or terms) of a
    !   complete polynomial of given degree and number of variables,
    !   evaluated at a given *single* point X in R^n.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), target :: x
        !*  Point X where basis functions should be evaluated
    integer, intent(in) :: k
        !*  (Maximum) degree of complete polynomial
    real (PREC), intent(out), dimension(:) :: basis
        !*  Basis functions evaluated at X
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    real (PREC), dimension(:,:), allocatable :: x2d, basis2d

    allocate (x2d(size(x),1), basis2d(size(basis),1))

    x2d(:,1) = x

    call polybasis_complete (x2d, k, basis2d, lstatus)

    if (lstatus == NF_STATUS_OK) then
        basis(:) = basis2d(:,1)
    end if
    if (present(status)) status = lstatus

    deallocate (x2d, basis2d)

end subroutine



pure subroutine __APPEND(polybasis_jac,__PREC) (x, k, jac, status)
    !*  POLYBASIS_JAC_COMPLETE returns the Jacobian of the basis functions
    !   (or terms) of a complete polynomial of given degree and number of variables,
    !   evaluated at a given set of x.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: x
        !*  Point where Jacobian of basis functions should be evaluated
    integer, intent(in) :: k
        !*  (Maximum) degree of complete polynomial
    real (PREC), intent(out), dimension(:,:) :: jac
        !*  Basis function Jacobian evaluated at X, shaped m-by-n
        !   where m is the number of basis functins and n is the dimension
        !   of X.
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    integer, dimension(:,:), allocatable :: exponents
    integer :: nbasis, nd, i, j, ib
    real (PREC) :: bi, xp

    lstatus = NF_STATUS_OK

    call check_nonneg (k, k, "k", lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! number of points and dimensions
    nd = size(x)
    ! Number of basis functions
    nbasis = poly_complete_get_nterms (nd, k)

    if (size(jac,1) < nbasis .or. size(jac,2) < nd) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Array that contains all permutations of exponents that occur in
    ! given complete polynomial (one per column)
    allocate (exponents(nd, nbasis))

    ! obtain sorted permutations of exponents that sum to k=2,...,k
    call polyexponents_complete (nd, k, exponents, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    do j = 1, nd
        ! Derivative of constant term is zero for all dimensions
        jac(1,:) = 0.0
        do ib = 2, nbasis

            ! Check if the ib-th basis function contains variable x_j; if
            ! not, the derivative is zero and we cycle to the next term.
            if (exponents(j,ib) == 0) then
                jac(ib,j) = 0.0
                cycle
            end if

            ! ib-th basis function contains x_j with positive exponent
            bi = 1.0
            do i = 1, nd
                xp = exponents(i,ib)
                if (i == j) then
                    bi = bi * xp * x(i) ** (xp - 1.0_PREC)
                else
                    bi = bi * x(i) ** xp
                end if
            end do
            jac(ib, j) = bi
        end do
    end do

100 continue

    if (allocated(exponents)) deallocate (exponents)

    if (present(status)) status = lstatus

end subroutine
