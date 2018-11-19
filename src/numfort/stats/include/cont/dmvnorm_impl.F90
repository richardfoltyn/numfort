

subroutine __APPEND(rvs_check_input,__PREC) (obj, x, mean, cov, dim, tol, status)
    integer, parameter :: PREC = __PREC
    type (__APPEND(dmvnorm,__PREC)), intent(in) :: obj
    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:), contiguous, target, optional :: mean
    real (PREC), intent(in), dimension(:,:), contiguous, target, optional :: cov
    integer, intent(in), optional :: dim
    real (PREC), intent(in), optional :: tol
    type (status_t), intent(out) :: status

    integer :: ldim
    real (PREC), dimension(:), pointer, contiguous :: ptr_mean
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_cov

    call get_dist_params (obj, mean, cov, ptr_mean, ptr_cov)

    if (.not. associated(ptr_mean) .or. .not. associated(ptr_cov)) goto 100

    if (size(ptr_mean) < 1) goto 100
    if (size(ptr_mean) /= size(ptr_cov,1)) goto 100
    if (size(ptr_cov,1) /= size(ptr_cov,2)) goto 100

    ldim = 1
    if (present(dim)) ldim = dim
    if (ldim /= 1 .or. ldim /= 2) goto 100
    if (size(x, ldim) /= size(mean)) goto 100

    call check_positive (0.0_PREC, tol, 'tol', status)
    if (status /= NF_STATUS_OK) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG
end subroutine



subroutine __APPEND(get_dist_params,__PREC) (obj, mean, cov, ptr_mean, ptr_cov)
    integer, parameter :: PREC = __PREC
    type (__APPEND(dmvnorm,__PREC)), intent(in), target :: obj
    real (PREC), intent(in), dimension(:), optional, target, contiguous :: mean
    real (PREC), intent(in), dimension(:,:), optional, target, contiguous :: cov
    real (PREC), intent(inout), dimension(:), pointer, contiguous :: ptr_mean
    real (PREC), intent(inout), dimension(:,:), pointer, contiguous :: ptr_cov

    nullify (ptr_mean, ptr_cov)

    if (present(mean)) then
        ptr_mean => mean
    else if (allocated(obj%mean)) then
        ptr_mean => mean
    end if

    if (present(cov)) then
        ptr_cov => cov
    else if (allocated(obj%cov)) then
        ptr_cov => cov
    end if

end subroutine



subroutine __APPEND(dmvnorm_rvs_2d,__PREC) (obj, x, mean, cov, dim, &
        check_cov, tol, status)
    integer, parameter :: PREC = __PREC
    type (__APPEND(dmvnorm,__PREC)), intent(in) :: obj
    real (PREC), intent(out), dimension(:,:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous, optional :: mean
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: cov
    integer, intent(in), optional :: dim
    logical, intent(in), optional :: check_cov
    real (PREC), intent(in), optional :: tol
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    logical :: lcheck_cov
    integer :: ldim
    real (PREC), dimension(:), contiguous, pointer :: ptr_mean
    real (PREC), dimension(:,:), contiguous, pointer :: ptr_cov

    integer :: i, j, ndim, nobs
    real (PREC), dimension(:,:), allocatable :: rwork2d
    real (PREC) :: diff, ltol

    ! Arguments used for call to GESDD
    real (PREC), dimension(:), allocatable :: s
    real (PREC), dimension(:,:), allocatable :: vt
    real (PREC), dimension(1,1) :: u
    real (PREC), dimension(:), allocatable :: rwork
    integer, dimension(:), allocatable :: iwork
    integer :: info, lwork

    lstatus = NF_STATUS_OK
    call rvs_check_input (obj, x, mean, cov, dim, tol, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ldim = 1
    lcheck_cov = .true.
    ltol = 1.0e-8_PREC
    if (present(dim)) ldim = dim
    if (present(check_cov)) lcheck_cov = check_cov
    if (present(tol)) ltol = tol

    call get_dist_params (obj, mean, cov, ptr_mean, ptr_cov)

    ! dimension of mv-normally distributed random vector
    ndim = size(ptr_mean)

    allocate (s(ndim), vt(ndim,ndim))
    allocate (rwork(1))
    allocate (iwork(8*ndim))

    ! Workspace query for SVD decomposition
    call LAPACK_GESDD (jobz='O', m=ndim, n=ndim, a=ptr_cov, lda=ndim, s=s, &
        u=u, ldu=1, vt=vt, ldvt=ndim, work=rwork, lwork=-1, iwork=iwork, &
        info=info)

    lwork = int(rwork(1))
    deallocate (rwork)
    allocate (rwork(lwork))

    ! Compute SVD decomposition (only S, V)
    call LAPACK_GESDD (jobz='O', m=ndim, n=ndim, a=ptr_cov, lda=ndim, s=s, &
        u=u, ldu=1, vt=vt, ldvt=ndim, work=rwork, lwork=lwork, iwork=iwork, &
        info=info)

    if (info /= 0) then
        if (info < 0 ) then
            ! This should not happend since we checked input arguments above
            lstatus = NF_STATUS_INVALID_ARG
        else
            lstatus = NF_STATUS_NOT_CONVERGED
        end if
        goto 100
    end if

    ! Convert V' to V * diag(sqrt(s))
    ! We keep calling it VT, though...
    allocate (rwork2d(ndim,ndim), source=vt)
    do i = 1, ndim
        vt(:,i) = rwork2d(i,:) * sqrt(s(i))
    end do

    ! Check whether COV was a positive semi-definite matrix, if applicable
    if (lcheck_cov) then
        call BLAS_GEMM (transa='T', transb='N', m=ndim, n=ndim, k=ndim, &
            alpha=1.0_PREC, a=vt, lda=ndim, b=vt, ldb=ndim, &
            beta=0.0_PREC, c=rwork2d, ldc=ndim)

        diff = maxval(abs(cov - rwork2d))
        if (diff > ltol) then
            lstatus = NF_STATUS_INVALID_STATE
            goto 100
        end if
    end if

    deallocate (rwork2d)

    ! Draw independent standard-normal realizations
    allocate (rwork2d(ndim, nobs))
    do j = 1, size(x, 2)
        do i = 1, size(x, 1)
            rwork2d(i,j) = __APPEND(random_normal,__PREC) ()
        end do
    end do

    ! Transform to potentially correlated MV-normal variables by computing
    !   X = V * Z
    ! where Z is the array of iid. standard-normal variables.
    nobs = size(x, 3-ldim)

    if (ldim == 1) then
        ! Each column contains one realized vector
        call BLAS_GEMM (transa='N', transb='N', m=ndim, n=nobs, k=ndim, &
            alpha=1.0_PREC, a=vt, lda=ndim, b=rwork2d, ldb=ndim, &
            beta=0.0_PREC, c=x, ldc=ndim)
    else
        ! Each row contains one realized vector, compute
        !   X' = Z' * V'
        ! instead.
        call BLAS_GEMM (transa='T', transb='T', m=nobs, n=ndim, k=ndim, &
            alpha=1.0_PREC, a=rwork2d, lda=ndim, b=vt, ldb=ndim, &
            beta=0.0_PREC, c=x, ldc=nobs)
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine
