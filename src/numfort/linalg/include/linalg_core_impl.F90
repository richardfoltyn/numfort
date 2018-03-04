

subroutine __APPEND(inv,__PREC) (A, Ainv, work, rwork, iwork, status)
    !*  INV computes the inverse of a matrix using the LAPACK routines
    !   GETRF and GETRI.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:), contiguous :: A
        !*  Input matrix to be inversted
    real (PREC), intent(out), dimension(:,:), contiguous :: Ainv
        !*  On exit, contains inverse of A if routine exits without errors.
    type (__APPEND(workspace,__PREC)), intent(in out), optional, target :: work
        !*  Workspace object (optional)
    real (PREC), intent(in out), dimension(:), optional, contiguous, target :: rwork
        !*  Working array for reals. If present, takes precedence over WORK.
    integer, intent(in out), dimension(:), optional, contiguous, target :: iwork
        !*  Working array for integers. If present, takes precedence over WORK.
    type (status_t), intent(out), optional :: status
        !*  Exit status (optoinal)

    type (__APPEND(workspace,__PREC)), pointer :: ptr_ws

    ! Arguments to LAPACK routines
    integer, dimension(:), pointer, contiguous :: ptr_ipiv
    real (PREC), dimension(:), pointer, contiguous :: ptr_work
    integer :: info
    integer :: n, lwork, lda, m
    real (PREC), dimension(1) :: rwork1
    integer, dimension(1) :: iwork1

    type (status_t) :: lstatus

    nullify (ptr_work, ptr_ipiv, ptr_ws)
    
    lstatus = NF_STATUS_OK

    if (.not. shape_equal (A, Ainv) .or. size(A,1) /= size(A,2)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    n = size(A,1)

    if (present(rwork)) then
        if (size(rwork) < n) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    if (present(iwork)) then
        if (size(iwork) < n) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A

    ! Workspace query for GETRI
    lda = n
    lwork = -1
    call LAPACK_GETRI (n, Ainv, lda, iwork1, rwork1, lwork, info)
    if (info /= 0) then
        lstatus = NF_STATUS_INVALID_STATE
        goto 100
    end if

    ! Optimal array size stored in first element of RWORK1
    lwork = int(rwork1(1))

    if (present(rwork) .and. present(iwork)) then
        if (lwork < size(rwork)) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
        ptr_work(1:lwork) => rwork
        ptr_ipiv(1:n) => iwork
    else
        call assert_alloc_ptr (work, ptr_ws)
        call workspace_reset (ptr_ws)

        call assert_alloc (ptr_ws, nrwrk=lwork, niwrk=n)

        call workspace_get_ptr (ptr_ws, n, ptr_work)
        call workspace_get_ptr (ptr_ws, n, ptr_ipiv)
    end if

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



subroutine __APPEND(det,__PREC) (A, d, work, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:), contiguous :: A
    real (PREC), intent(out) :: d
    type (__APPEND(workspace,__PREC)), intent(in out), optional, target :: work
    type (status_t), intent(out), optional :: status

    type (__APPEND(workspace,__PREC)), pointer :: ptr_ws
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
