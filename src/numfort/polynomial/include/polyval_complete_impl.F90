

pure subroutine __APPEND(polybasis_check_input_array,__PREC) (x, k, basis, &
        exponents, status)
    !*  POLYBASIS_CHECK_INPUT performs input checks for user-provided
    !   arguments for POLYBASIS_COMPLETE.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: x
    integer, intent(in) :: k
    real (PREC), intent(in), dimension(:,:) :: basis
    integer, intent(in), dimension(:,:), optional :: exponents
    type (status_t), intent(out), optional :: status

    integer :: npoints, nd, nbasis

    status = NF_STATUS_OK

    npoints = size(x, 2)
    nd = size(x, 1)

    call check_nonneg (k, k, "k", status=status)
    if (status /= NF_STATUS_OK) goto 100

    nbasis = poly_complete_get_nterms (nd, k)

    if (size(basis, 1) < nbasis .or. size(basis, 2) < npoints) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(exponents)) then
        if (size(exponents, 1) /= nd .or. size(exponents, 2) /= nbasis) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

100 continue

end subroutine



pure subroutine __APPEND(polybasis_check_input_scalar,__PREC) (x, k, basis, &
        exponents, status)
    !*  POLYBASIS_CHECK_INPUT performs input checks for user-provided
    !   arguments for POLYBASIS_COMPLETE.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    integer, intent(in) :: k
    real (PREC), intent(in), dimension(:) :: basis
    integer, intent(in), dimension(:,:), optional :: exponents
    type (status_t), intent(out), optional :: status

    integer :: nd, nbasis

    status = NF_STATUS_OK

    nd = size(x)

    call check_nonneg (k, k, "k", status=status)
    if (status /= NF_STATUS_OK) goto 100

    nbasis = poly_complete_get_nterms (nd, k)

    if (size(basis) < nbasis) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(exponents)) then
        if (size(exponents, 1) /= nd .or. size(exponents, 2) /= nbasis) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

100 continue

end subroutine


pure subroutine __APPEND(polybasis,__PREC) (x, k, basis, exponents, status)
    !*  POLYBASIS_COMPLETE returns the basis functions (or terms) of a
    !   complete polynomial of given degree and number of variables,
    !   evaluated at a given set of x.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  Set of x where basis functions should be evaluated
        !   (each column contains on point of dimension N)
    integer, intent(in) :: k
        !*  (Maximum) degree of complete polynomial
    real (PREC), intent(out), dimension(:,:), contiguous :: basis
        !*  Basis functions evaluated at X
    integer, intent(in), dimension(:,:), optional, contiguous :: exponents
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    integer, dimension(:,:), allocatable :: lexponents
    integer :: nbasis, nd, j, nx

    lstatus = NF_STATUS_OK

    call polybasis_check_input (x, k, basis, exponents, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! number of points and dimensions
    nd = size(x, 1)
    nx = size(x, 2)
    ! No points given, hence nothing to do
    if (nx == 0) goto 100

    ! Basis functions for 0-dimensional data: in this case we should just
    ! return the constant term and exit immediately.
    if (nd == 0) then
        basis = 0.0_PREC
        basis(1,:) = 1.0_PREC
        goto 100
    end if

    ! Number of basis functions
    nbasis = poly_complete_get_nterms (nd, k)

    if (.not. present(exponents)) then
        ! Array of exponents for each basis function term needs to be
        ! computed locally.

        ! Array that contains all permutations of exponents that occur in
        ! given complete polynomial (one per column).
        allocate (lexponents(nd, nbasis))

        ! obtain sorted permutations of exponents that sum to k=2,...,k
        call polyexponents_complete (nd, k, lexponents, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        do j = 1, nx
            call polybasis_complete_impl (x(:,j), lexponents, basis(:,j))
        end do
    else
        do j = 1, nx
            call polybasis_complete_impl (x(:,j), exponents, basis(:,j))
        end do
    end if


100 continue

    if (allocated(lexponents)) deallocate (lexponents)

    if (present(status)) status = lstatus

end subroutine


pure subroutine __APPEND(polybasis_scalar,__PREC) (x, k, basis, exponents, status)
    !*  POLYBASIS_COMPLETE_SCALAR returns the basis functions (or terms) of a
    !   complete polynomial of given degree and number of variables,
    !   evaluated at a given *single* point X in R^n.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Point X where basis functions should be evaluated
    integer, intent(in) :: k
        !*  (Maximum) degree of complete polynomial
    real (PREC), intent(out), dimension(:), contiguous :: basis
        !*  Basis functions evaluated at X
    integer, intent(in), dimension(:,:), optional, contiguous :: exponents
        !*  Array that contains exponents for each variable (ie. dimension)
        !   for every basis function term. Shaped [size(X), size(BASIS)]
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    integer, dimension(:,:), allocatable :: lexponents
    integer :: nbasis, nd

    lstatus = NF_STATUS_OK

    call polybasis_check_input (x, k, basis, exponents, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    nd = size(x)
    if (nd == 0.0) then
        basis = 0.0_PREC
        basis(1) = 1.0_PREC
        goto 100
    end if

    ! Number of basis functions
    nbasis = poly_complete_get_nterms (nd, k)

    if (.not. present(exponents)) then
        ! Array of exponents for each basis function term needs to be
        ! computed locally.

        ! Array that contains all permutations of exponents that occur in
        ! given complete polynomial (one per column).
        allocate (lexponents(nd, nbasis))

        ! obtain sorted permutations of exponents that sum to k=2,...,k
        call polyexponents_complete (nd, k, lexponents, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call polybasis_complete_impl (x, lexponents, basis)
    else
        call polybasis_complete_impl (x, exponents, basis)
    end if

100 continue

    if (allocated(lexponents)) deallocate (lexponents)
    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(polybasis_impl,__PREC) (x, exponents, basis)
    !*  POLYBASIS_COMPLETE returns the basis functions (or terms) of a
    !   complete polynomial of given degree and number of variables,
    !   evaluated at a given set of x.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Set of x where basis functions should be evaluated
        !   (each column contains on point of dimension N)
    integer, intent(in), dimension(:,:), optional, contiguous :: exponents
        !*  Array that contains exponents for each variable (ie. dimension)
        !   for every basis function term. Shaped [size(X), size(BASIS)]
    real (PREC), intent(out), dimension(:), contiguous :: basis
        !*  Basis functions evaluated at X

    integer :: nd, nbasis, i, ib
    real (PREC) :: bi

    ! number of points and dimensions
    nd = size(x)
    nbasis = size(exponents, 2)

    do ib = 1, nbasis
        bi = 1.0
        do i = 1, nd
            bi = bi * x(i) ** exponents(i, ib)
        end do
        basis(ib) = bi
    end do

end subroutine



pure subroutine __APPEND(polybasis_jac,__PREC) (x, k, jac, exponents, status)
    !*  POLYBASIS_JAC_COMPLETE returns the Jacobian of the basis functions
    !   (or terms) of a complete polynomial of given degree and number of variables,
    !   evaluated at a given set of x.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Point where Jacobian of basis functions should be evaluated
    integer, intent(in) :: k
        !*  (Maximum) degree of complete polynomial
    real (PREC), intent(out), dimension(:,:), contiguous :: jac
        !*  Basis function Jacobian evaluated at X, shaped m-by-n
        !   where m is the number of basis functins and n is the dimension
        !   of X.
    integer, intent(in), dimension(:,:), optional, contiguous :: exponents
        !*  Array that contains exponents for each variable (ie. dimension)
        !   for every basis function term. Shaped [size(X), NBASIS]
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    integer, dimension(:,:), allocatable :: lexponents
    integer :: nbasis, nd

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

    if (present(exponents)) then
        if (size(exponents,1) /= nd .or. size(exponents,2) /= nbasis) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    if (present(exponents)) then
        call polybasis_jac_complete_impl (x, exponents, jac)
    else
         ! Array that contains all permutations of exponents that occur in
        ! given complete polynomial (one per column)
        allocate (lexponents(nd, nbasis))

        ! obtain sorted permutations of exponents that sum to k=2,...,k
        call polyexponents_complete (nd, k, lexponents, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call polybasis_jac_complete_impl (x, lexponents, jac)
    end if

100 continue

    if (allocated(lexponents)) deallocate (lexponents)
    if (present(status)) status = lstatus

end subroutine


pure subroutine __APPEND(polybasis_jac_impl,__PREC) (x, exponents, jac)
    !*  POLYBASIS_JAC_COMPLETE returns the Jacobian of the basis functions
    !   (or terms) of a complete polynomial of given degree and number of variables,
    !   evaluated at a given set of x.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Point where Jacobian of basis functions should be evaluated
    real (PREC), intent(out), dimension(:,:), contiguous :: jac
        !*  Basis function Jacobian evaluated at X, shaped m-by-n
        !   where m is the number of basis functins and n is the dimension
        !   of X.
    integer, intent(in), dimension(:,:), contiguous :: exponents
        !*  Array that contains exponents for each variable (ie. dimension)
        !   for every basis function term. Shaped [size(X), NBASIS]

    integer :: nbasis, nd, i, j, ib
    real (PREC) :: bi, xp

    nd = size(x)
    nbasis = size(exponents, 2)

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

end subroutine

