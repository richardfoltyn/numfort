

subroutine inv_work_query_matrix (A, lwork, liwork, status)
    !*  INV_WORK_QUERY returns the optimal workspace sizes for real and
    !   integer working arrays needed by the INV routine.
    real (PREC), intent(in), dimension(:,:) :: A
        !*  Matrix that should be inverted.
    integer, intent(out) :: lwork
        !*  Optimal REAL working array size
    integer, intent(out) :: liwork
        !*  INTEGER working array size
    type (status_t), intent(out), optional :: status
        !*  Exit status.

    call inv_work_query (0.0_PREC, size(A,1), lwork, liwork, status)

end subroutine



subroutine inv_work_query_dims (dummy, n, lwork, liwork, status)
    !*  INV_WORK_QUERY returns the optimal workspace sizes for real and
    !   integer working arrays needed by the INV routine.
    real (PREC), intent(in) :: dummy
    integer, intent(in) :: n
    !*  Number of rows (columns) of square matrix to be inverted.
    integer, intent(out) :: lwork
    !*  Optimal REAL working array size
    integer, intent(out) :: liwork
    !*  INTEGER working array size
    type (status_t), intent(out), optional :: status
    !*  Exit status.

    type (status_t) :: lstatus
    integer :: lda, info
    integer, dimension(1) :: idummy1d
    real (PREC), dimension(1) :: rdummy1d
    real (PREC), dimension(1,1) :: rdummy2d

    lstatus = NF_STATUS_OK

    liwork = -1

    ! Workspace query for GETRI
    lda = n
    lwork = -1
    call LAPACK_GETRI (n, rdummy2d, lda, idummy1d, rdummy1d, lwork, info)
    if (info /= 0) then
        lstatus = NF_STATUS_INVALID_STATE
        goto 100
    end if

    ! Optimal array size stored in first element of RWORK1
    lwork = int(rdummy1d(1))
    liwork = n

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine inv (A, Ainv, work, rwork, iwork, status)
    !*  INV computes the inverse of a matrix using the LAPACK routines
    !   GETRF and GETRI.
    real (PREC), intent(in), dimension(:,:), contiguous :: A
        !*  Input matrix to be inversted
    real (PREC), intent(out), dimension(:,:), contiguous :: Ainv
        !*  On exit, contains inverse of A if routine exits without errors.
    type (workspace), intent(inout), optional, target :: work
        !*  Workspace object (optional)
    real (PREC), intent(inout), dimension(:), optional, contiguous, target :: rwork
        !*  Working array for reals. If present, takes precedence over WORK.
    integer, intent(inout), dimension(:), optional, contiguous, target :: iwork
        !*  Working array for integers. If present, takes precedence over WORK.
    type (status_t), intent(out), optional :: status
        !*  Exit status (optoinal)

    type (workspace), pointer :: ptr_ws

    ! Arguments to LAPACK routines
    integer, dimension(:), pointer, contiguous :: ptr_ipiv
    real (PREC), dimension(:), pointer, contiguous :: ptr_work
    integer :: info
    integer :: n, lwork, lda, m, liwork

    type (status_t) :: lstatus

    nullify (ptr_work, ptr_ipiv, ptr_ws)
    
    lstatus = NF_STATUS_OK

    if (.not. shape_equal (A, Ainv) .or. size(A,1) /= size(A,2)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    n = size(A,1)

    if (present(rwork) .and. present(iwork)) then
        ! Enforce minimum size requirements
        if (n > size(rwork)) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if

        if (n > size(iwork)) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if

        lwork = size(rwork)
        ptr_work(1:lwork) => rwork
        ptr_ipiv(1:n) => iwork(1:n)
    else
        call assert_alloc_ptr (work, ptr_ws)
        call workspace_reset (ptr_ws)
        
        call inv_work_query (A, lwork, liwork)

        call assert_alloc (ptr_ws, nrwrk=lwork, niwrk=liwork)

        call workspace_get_ptr (ptr_ws, lwork, ptr_work)
        call workspace_get_ptr (ptr_ws, liwork, ptr_ipiv)
    end if

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    lda = n

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    m = n

    call LAPACK_GETRF (m, n, Ainv, lda, ptr_ipiv, info)

    if (info < 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    else if (info > 0) then
        ! Factor U is exactly singular
        lstatus = NF_STATUS_INVALID_STATE
        goto 100
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call LAPACK_GETRI (n, Ainv, lda, ptr_ipiv, ptr_work, lwork, info)

    if (info < 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    else if (info > 0) then
        ! Factor U is exactly singular, inverse could not be computed
        lstatus = NF_STATUS_INVALID_STATE
        goto 100
    end if

100 continue

    if (present(status)) status = lstatus

    ! Deallocate WORKSPACE object if it was allocated locally
    call assert_dealloc_ptr (work, ptr_ws)

end subroutine



subroutine det (A, d, work, status)
    real (PREC), intent(in), dimension(:,:), contiguous :: A
    real (PREC), intent(out) :: d
    type (workspace), intent(inout), optional, target :: work
    type (status_t), intent(out), optional :: status

    type (workspace), pointer :: ptr_ws
    integer :: nrwrk, niwrk
    integer, dimension(2) :: shp2d
    integer :: pow, i

    ! Arguments to LAPACK routines
    integer, dimension(:), pointer, contiguous :: ptr_ipiv
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_a
    integer :: info
    integer :: n, m, lda

    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    d = 0.0_PREC

    m = size(A,1)
    n = size(A,2)
    nrwrk = m*n
    niwrk = min(m,n)

    call assert_alloc_ptr (work, ptr_ws)
    call workspace_reset (ptr_ws)

    call assert_alloc (ptr_ws, nrwrk=nrwrk, niwrk=niwrk)

    shp2d(1) = m
    shp2d(2) = n
    call workspace_get_ptr (ptr_ws, shp2d, ptr_a)
    call workspace_get_ptr (ptr_ws, niwrk, ptr_ipiv)

    lda = n
    m = n
    call LAPACK_GETRF (m, n, ptr_a, lda, ptr_ipiv, info)

    if (info < 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    else if (info > 0) then
        ! Factor U is exactly singular
        lstatus = NF_STATUS_INVALID_STATE
        goto 100
    end if

    pow = 0
    d = 1.0_PREC

    do i = 1, n
        d = d * ptr_a(i,i)
        if (ptr_ipiv(i) /= i) pow = pow + 1
    end do

    ! Correct sign
    d = d * (-1) ** pow

100 continue
    if (present(status)) status = lstatus

    ! Deallocate WORKSPACE object if it was allocated locally
    call assert_dealloc_ptr (work, ptr_ws)

end subroutine



subroutine gram (a, g, trans, status)
    !*  GRAM computes the Gram matrix G = A'A, or alternatively
    !   the left Gram matrix G = AA'.
    real (PREC), dimension(:,:), intent(in), contiguous :: a
        !*   Matrix A.
    real (PREC), dimension(:,:), intent(out), contiguous :: g
        !*  Gram matrix
    character (1), intent(in), optional :: trans
        !*  If present and equal to T, the Gram matrix G = A'A is returned.
        !   If present and equal to N, the left Gram matrix G = AA'
        !   is computed (default: trans = 'T').
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    type (status_t) :: lstatus
    integer :: k, n, lda, ldc, i, j
    character (1) :: ltrans
    character (*), parameter :: NAME = 'GRAM'
    character (*), parameter :: UPLO = 'U'
    real (PREC), parameter :: alpha = 1.0_PREC
    real (PREC), parameter :: beta = 0.0_PREC

    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    ltrans = 'T'
    if (present(trans)) ltrans = trans

    call lower (ltrans)

    call check_cond (ltrans == 't' .or. ltrans == 'n', NAME, &
        'TRANS: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Determine dimensions of the problem, using the same naming scheme as SYRK
    if (ltrans == 't') then
        k = size(a, 1)
        n = size(a, 2)
    else
        k = size(a, 2)
        n = size(a, 1)
    end if

    call check_cond (size(g,1) == n .and. size(g,2) == n, NAME, &
        'G: Non-conformable array size', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Implementation ---

    lda = size(a, 1)
    ldc = n

    call BLAS_SYRK (uplo, ltrans, n, k, alpha, a, lda, beta, g, ldc)

    ! Reflect matrix G along its diagonal. SYRK writes data to upper triangular
    ! portion of G.

    do j = 1, n
        do i = j+1, n
            g(i,j) = g(j,i)
        end do
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine
