
! Define LAPACK routines for given precision

#if __PREC == real64
#define __GESVD DGESVD
#define __GEMM DGEMM
#define __GELSD DGELSD
#define __GESDD DGESDD
#elif __PREC == real32
#define __GESVD SGESVD
#define __GEMM SGEMM
#define __GELSD SGELSD
#define __GESDD SGESDD
#endif

!-------------------------------------------------------------------------------
! OLS

pure subroutine __APPEND(ols_check_input,__PREC) (x, y, beta, add_const, trans_x, &
        rcond, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    real (PREC), intent(in), dimension(:,:) :: beta
    logical, intent(in) :: add_const
    logical, intent(in) :: trans_x
    real (PREC), intent(in) :: rcond
    type (status_t), intent(out) :: status

    integer :: nvars, nobs, nconst, nlhs

    status = NF_STATUS_INVALID_ARG

    nconst = 0
    if (add_const) nconst = 1

    if (.not. trans_x) then
        nobs = size(x, 1)
        nvars = size(x, 2)
    else
        nobs = size(x, 2)
        nvars = size(x, 1)
    end if

    nlhs = size(y, 2)

    if (nobs < nvars + nconst) return
    if (size(y, 1) /= nobs) return
    if (size(beta, 2) /= nlhs) return
    if (size(beta, 1) /= (nvars + nconst)) return
    if (rcond < 0) return

    status = NF_STATUS_OK

end subroutine


subroutine __APPEND(ols_2d,__PREC) (x, y, beta, add_const, trans_x, rcond, &
        rank, status)

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    real (PREC), intent(out), dimension(:,:) :: beta
    logical, intent(in), optional :: add_const
    logical, intent(in), optional :: trans_x
    real (PREC), intent(in), optional :: rcond
    integer, intent(out), optional :: rank
    type (status_t), intent(out), optional :: status

    logical :: ladd_const, ltrans_x
    real (PREC) :: lrcond
    integer :: nobs, nvars, ncoefs, nconst, nlhs, i
    type (status_t) :: lstatus

    ! GELSD arguments
    real (PREC), dimension(:), allocatable :: work, sval
    integer, dimension(:), allocatable :: iwork
    real (PREC), dimension(:,:), allocatable :: lx, lhs
    integer :: lrank
    integer :: lwork, info, m, n, nrhs, lda, ldb, mn, liwork

    ladd_const = .true.
    if (present(add_const)) ladd_const = add_const

    ltrans_x = .false.
    if (present(trans_x)) ltrans_x = trans_x

    ! Use default value from LAPACK95 GELS wrapper
    lrcond = 100 * epsilon(1.0_PREC)
    if (present(rcond)) lrcond = rcond

    call ols_check_input (x, y, beta, ladd_const, ltrans_x, lrcond, lstatus)

    nvars = size(x, 2)
    if (ltrans_x) nvars = size(x, 1)
    nobs = size(y, 1)
    nlhs = size(y, 2)
    nconst = 0
    if (ladd_const) nconst = 1
    ncoefs = nvars + nconst

    ! Allocate (possibly transposed) array to store X variables; will be
    ! overwritten by GELSD, so we have to allocate in any case.
    ! Add constant as requested.
    allocate (lx(nobs,ncoefs))
    if (ltrans_x) then
        forall (i=1:nvars) lx(:,nconst+i) = x(i,:)
    else
        lx(:, 1+nconst:ncoefs) = x
    end if
    if (ladd_const) lx(:, 1) = 1.0_PREC

    ! Set up GELSD call
    m = nobs
    n = ncoefs
    mn = max(1, min(m, n))
    lda = nobs
    ldb = nobs
    ! GELSD solves Ax=b, we solve y = Xb, so our nlhs = nrhs in GELSD
    nrhs = nlhs
    allocate (sval(mn))
    ! Dummy argument b will be overwritten, create copy
    allocate (lhs, source=y)

    ! 1. Workspace query
    lwork = -1
    allocate (work(1))
    allocate (iwork(1))

    call __GELSD (m, n, nrhs, lx, lda, lhs, ldb, sval, lrcond, lrank, work, &
        lwork, iwork, info)

    if (info < 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! 2. Allocate optimal / minimal work arrays
    lwork = int(work(1))
    liwork = iwork(1)
    
    deallocate (work, iwork)
    allocate (work(lwork), iwork(liwork))

    ! 3. Actual call to perform LS
    call __GELSD (m, n, nrhs, lx, lda, lhs, ldb, sval, lrcond, lrank, work, &
        lwork, iwork, info)

    ! Copy coefficients
    forall (i=1:nlhs) beta(1:ncoefs,i) = lhs(1:ncoefs,i)
    ! Copy over optional output arguments
    if (present(rank)) rank = lrank

100 continue
    if (present(status)) status = lstatus

end subroutine

subroutine __APPEND(ols_1d,__PREC) (x, y, beta, add_const, trans_x, rcond, &
        rank, status)

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:), target :: y
    real (PREC), intent(out), dimension(:), target :: beta
    logical, intent(in), optional :: add_const
    logical, intent(in), optional :: trans_x
    real (PREC), intent(in), optional :: rcond
    integer, intent(out), optional :: rank
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer :: ptr_y, ptr_beta
    integer :: nobs, ncoefs

    nobs = size(y)
    ncoefs = size(beta)

    ptr_y(1:nobs,1:1) => y
    ptr_beta(1:ncoefs,1:1) => beta

    call ols (x, ptr_y, ptr_beta, add_const, trans_x, rcond, rank, status)

end subroutine

!-------------------------------------------------------------------------------
! PCA (Principal component analysis)

pure subroutine __APPEND(pca_check_input,__PREC) (x, scores, ncomp, s, &
        loadings, mean_x, std_x, propvar, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:,:) :: scores
    integer, intent(in) :: ncomp
    real (PREC), intent(in), dimension(:), optional :: s
    real (PREC), intent(in), dimension(:,:), optional :: loadings
    real (PREC), intent(in), dimension(:), optional :: mean_x, std_x
    real (PREC), intent(in), dimension(:), optional :: propvar
    type (status_t), intent(out) :: status

    integer :: nobs, nvars

    nobs = size(x, 1)
    nvars = size(x, 2)

    status = NF_STATUS_INVALID_ARG

    if (nobs < 1 .or. nobs < nvars) return
    if (ncomp < 1 .or. ncomp > nvars) return
    if (size(scores,1) < nobs .or. size(scores,2) < ncomp) return

    if (present(s)) then
        if (size(s) < ncomp) return
    end if

    if (present(loadings)) then
        if (size(loadings, 1) < nvars .or. size(loadings,2) < ncomp) return
    end if

    if (present(mean_x)) then
        if (size(mean_x) < nvars) return
    end if

    if (present(std_x)) then
        if (size(std_x) < nvars) return
    end if

    ! Proportional of variance explained by each of the first ncomp PCs
    if (present(propvar)) then
        if (size(propvar) < ncomp) return
    end if

    status = NF_STATUS_OK

end subroutine


subroutine __APPEND(pca,__PREC) (x, scores, ncomp, center, scale, &
        sval, loadings, mean_x, std_x, propvar, status)

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(out), dimension(:,:), contiguous :: scores
    integer, intent(in) :: ncomp
    logical, intent(in), optional :: center
    logical, intent(in), optional :: scale
    real (PREC), intent(out), dimension(:), optional :: sval
    real (PREC), intent(out), dimension(:,:), optional :: loadings
    real (PREC), intent(out), dimension(:), optional :: mean_x
    real (PREC), intent(out), dimension(:), optional :: std_x
    real (PREC), intent(out), dimension(:), optional :: propvar
    type (status_t), intent(out), optional :: status

    logical :: lscale, lcenter
    type (status_t) :: lstatus
    real (PREC), dimension(:,:), allocatable :: x_n, x_tmp
    real (PREC), dimension(:), allocatable :: work
    integer :: nvars, nobs, i
    real (PREC), dimension(:), allocatable :: lmean_x, lstd_x

    ! Argument for GESDD
    real (PREC), dimension(:), allocatable :: ls
    real (PREC), dimension(:,:), allocatable :: vt, u
    real (PREC), dimension(1) :: qwork
    integer, dimension(:), allocatable :: iwork
    character (1), parameter :: jobz = 'A'
    integer :: lda, ldu, ldvt, lwork, info, m, n, mn

    ! Arguments for GEMM
    real (PREC), dimension(:,:), allocatable :: vk
    character (1) :: transa, transb
    real (PREC), parameter :: alpha = 1.0_PREC, beta = 0.0_PREC
    integer :: ldb, ldc, k

    call pca_check_input (x, scores, ncomp, sval, loadings, mean_x, std_x, propvar, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    lcenter = .true.
    lscale = .true.
    if (present(center)) lcenter = center
    if (present(scale)) lscale = scale

    nobs = size(x, 1)
    nvars = size(x, 2)

    allocate (x_n, source=x)

    allocate (lmean_x(nvars), lstd_x(nvars))

    ! Normalize data in place
    if (lcenter .or. lscale) then
        call normalize (x_n, m=lmean_x, s=lstd_x, scale=lscale, center=lcenter, dof=0, dim=1)
    end if

    m = nobs
    n = nvars
    lda = nobs
    ldvt = n
    ldu = m
    mn = max(1, min(m, n))

    allocate (iwork(8*mn))
    allocate (vt(n, n))
    allocate (ls(mn))
    allocate (u(m,m))
    ! Contents of x_n would be overwritten by GESDD, make copy
    allocate (x_tmp, source=x_n)

    ! workspace query
    lwork = -1
    call __GESDD (jobz, m, n, x_tmp, lda, ls, u, ldu, vt, ldvt, qwork, lwork, iwork, info)

    ! Recover minimal work space size
    if (info /= 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if
    lwork = int(qwork(1))

    ! perform actual SVD
    allocate (work(lwork))

    call __GESDD (jobz, m, n, x_tmp, lda, ls, u, ldu, vt, ldvt, work, lwork, iwork, info)
    deallocate (x_tmp)

    if (info /= 0) then
        if (info < 0) then
            ! This should not happen since we checked arguments before
            lstatus = NF_STATUS_INVALID_ARG
        else
            lstatus = NF_STATUS_NOT_CONVERGED
        end if
        goto 100
    end if

    ! Copy over first ncomp components and transpose
    allocate (vk(nvars, ncomp))
    forall (i=1:ncomp) vk(:,i) = vt(i, :)

    ! Compute the first k principal components (scores) as PC_k = X V_k
    m = nobs
    n = ncomp
    k = nvars
    lda = nobs
    ldb = nvars
    ldc = nobs
    transa = 'N'
    transb = 'N'
    call __GEMM (transa, transb, m, n, k, alpha, x_n, lda, vk, ldb, beta, scores, ldc)

    if (present(sval)) sval(1:ncomp) = ls(1:ncomp)
    if (present(loadings)) then
        forall (i=1:ncomp) loadings(1:nvars,i) = vk(:, i)
    end if

    ! compute vector of variance share explained by each PC.
    ! Variance share is given by relative size of each eigenvalue (ie. squared
    ! singular value) relative to sum of all EVs.
    if (present(propvar)) then
        propvar(1:ncomp) = ls(1:ncomp) ** 2 / sum(ls ** 2)
    end if

    if (present(mean_x) .and. lcenter) mean_x = lmean_x
    if (present(std_x) .and. lscale) std_x = lstd_x

100 continue
    if (present(status)) status = lstatus

end subroutine

! ------------------------------------------------------------------------------
! PCR (Principal component regression)


pure subroutine __APPEND(pcr_check_input,__PREC) (lhs, scores, sval, loadings, &
        coefs, mean_x, std_x, add_const, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: lhs
    real (PREC), intent(in), dimension(:,:) :: scores
    real (PREC), intent(in), dimension(:) :: sval
    real (PREC), intent(in), dimension(:,:) :: loadings
    real (PREC), intent(in), dimension(:,:) :: coefs
    real (PREC), intent(in), dimension(:), optional :: mean_x
    real (PREC), intent(in), dimension(:), optional :: std_x
    logical, intent(in) :: add_const
    type (status_t), intent(out) :: status

    integer :: nobs, nvars, ncomp, nlhs, nconst

    nconst = 0
    if (add_const) nconst = 1

    nobs = size(lhs, 1)
    nlhs = size(lhs, 2)
    ncomp = size(scores, 2)
    nvars = size(coefs, 1) - nconst

    status = NF_STATUS_INVALID_ARG

    if (nobs < 1 .or. nobs < ncomp) return
    if (size(scores, 1) /= nobs) return
    if (size(sval) < ncomp) return
    if (any(shape(loadings) /= [nvars, ncomp])) return
    if (any(shape(coefs) /= [nvars + nconst, nlhs])) return
    if (present(mean_x)) then
        if (size(mean_x) < nvars) return
    end if
    if (present(std_x)) then
        if (size(std_x) < nvars) return
    end if

    status = NF_STATUS_OK

end subroutine

subroutine __APPEND(pcr_2d,__PREC) (lhs, scores, sval, loadings, coefs, mean_x, std_x, &
        center, add_const, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: lhs
    real (PREC), intent(in), dimension(:,:) :: scores
    real (PREC), intent(in), dimension(:) :: sval
    real (PREC), intent(in), dimension(:,:) :: loadings
    real (PREC), intent(out), dimension(:,:) :: coefs
    real (PREC), intent(in), dimension(:), optional :: mean_x
    real (PREC), intent(in), dimension(:), optional :: std_x
    logical, intent(in), optional :: add_const
    logical, intent(in), optional :: center
    type (status_t), intent(out), optional :: status

    logical :: ladd_const, lcenter
    integer :: nvars, nobs, ncomp, nlhs
    real (PREC), dimension(:), pointer :: ptr_mean_x, ptr_std_x

    type (status_t) :: lstatus
    integer :: i, nconst

    real (PREC), dimension(:,:), allocatable :: y_n, PCy, lcoefs
    real (PREC), dimension(:), allocatable :: mean_y
    ! Arguments to GEMM
    character (1) :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real (PREC), parameter :: alpha = 1.0_PREC, beta = 0.0_PREC

    ladd_const = .true.
    lcenter = .true.
    if (present(add_const)) ladd_const = add_const
    if (present(center)) lcenter = center

    call pcr_check_input (lhs, scores, sval, loadings, coefs, mean_x, std_x, &
        ladd_const, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    nconst = 0
    if (ladd_const) nconst = 1
    nvars = size(coefs, 1) - nconst
    ncomp = size(loadings, 2)
    nobs = size(lhs, 1)
    nlhs = size(lhs, 2)

    ! sample mean and standard deviation of original x variables
    call assert_alloc_ptr (mean_x, nvars, ptr_mean_x)
    call assert_alloc_ptr (std_x, nvars, ptr_std_x)

    if (.not. present(mean_x)) ptr_mean_x = 0.0_PREC
    if (.not. present(std_x)) ptr_std_x = 1.0_PREC

    ! copy over dependent variable, will be overwritten by GELS
    allocate (y_n(nobs, nlhs), source=lhs)
    allocate (mean_y(nlhs), source=0.0_PREC)
    if (lcenter) then
        call normalize (y_n, m=mean_y, dim=1, center=.true., scale=.false.)
    end if

    ! Compute regression coefficients of regression y_n on PC
    ! 1. Compute PC'y
    m = ncomp
    n = nlhs
    k = nobs
    lda = nobs
    ldb = nobs
    ldc = ncomp
    transa = 'T'
    transb = 'N'
    allocate (PCy(ncomp, nlhs))
    call __GEMM (transa, transb, m, n, k, alpha, scores, lda, y_n, ldb, beta, PCy, ldc)

    ! 2. Regression coefficients are given by (diag(sval**2))^{-1} PC'y,
    ! can be computed directly without solving equation system.
    do i = 1, nlhs
        PCy(1:ncomp, i) = PCy(1:ncomp, i) * (sval(1:ncomp) ** (-2.0_PREC))
    end do

    ! PCR coefficients derived from regressing on k PC are given by
    ! beta_pcr = V_k * beta_ols
    m = nvars
    n = nlhs
    k = ncomp
    lda = nvars
    ldb = ncomp
    ldc = nvars
    transa = 'N'
    transb = 'N'
    allocate (lcoefs(nvars, nlhs))
    call __GEMM (transa, transb, m, n, k, alpha, loadings, lda, PCy, ldb, beta, lcoefs, ldc)

    ! Scale back betas if explanatory variables were normalized
    do i = 1, nlhs
        lcoefs(:, i) = lcoefs(:, i) / std_x
    end do

    ! Add intercept if requested
    if (ladd_const) then
        do i = 1, nlhs
            ! set intercept such that sample means lie on regression line
            coefs(1, i) = mean_y(i) - dot_product(mean_x, lcoefs(:, i))
            coefs(2:nvars+1,i) = lcoefs(:, i)
        end do
    else
        forall (i=1:nlhs) coefs(1:nvars,i) = lcoefs(:,i)
    end if

100 continue
    call assert_dealloc_ptr (mean_x, ptr_mean_x)
    call assert_dealloc_ptr (std_x, ptr_std_x)
    if (present(status)) lstatus = status

end subroutine

subroutine __APPEND(pcr_1d,__PREC) (lhs, scores, sval, loadings, coefs, &
        mean_x, std_x, add_const, center, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), target :: lhs
    real (PREC), intent(in), dimension(:,:) :: scores
    real (PREC), intent(in), dimensioN(:) :: sval
    real (PREC), intent(in), dimension(:,:) :: loadings
    real (PREC), intent(out), dimension(:), target :: coefs
    real (PREC), intent(in), dimension(:), optional :: mean_x
    real (PREC), intent(in), dimension(:), optional :: std_x
    logical, intent(in), optional :: add_const
    logical, intent(in), optional :: center
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer :: ptr_lhs, ptr_coefs

    integer :: ncoefs, nobs

    nobs = size(lhs, 1)
    ncoefs = size(coefs, 1)

    ptr_lhs(1:nobs,1:1) => lhs
    ptr_coefs(1:ncoefs,1:1) => coefs

    call pcr (ptr_lhs, scores, sval, loadings, ptr_coefs, &
        mean_x, std_x, add_const, center, status)

end subroutine

! pure subroutine get_pcr_work_size (m, n, nrwork, niwork)
!     integer, intent(in) :: m, n
!     integer, intent(out) :: nrwork, niwork
!
!     integer :: smlsiz, nlvl, mn
!
!     mn = min(m, n)
!
!     ! Copied from Netlib's LAPACK95 wrappers
!     smlsiz = ilaenv (9, 'DGELSD', ' ', 0, 0, 0, 0)
!     nlvl = int(log(max(1, mn) / (smlsiz + 1.0d0)) / log(2.0d0)) + 1
!
!     nrwork = 12*mn + 2*mn*smlsiz + 8*mn*nlvl + mn*rhs + (smlsiz+1) ** 2
!
! end subroutine


#ifdef __GESVD
#undef __GESVD
#endif

#ifdef __GEMM
#undef __GEMM
#endif

#ifdef __GELSD
#undef __GELSD
#endif
