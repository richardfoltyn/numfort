
! Define LAPACK routines for given precision

!-------------------------------------------------------------------------------
! OLS

pure subroutine __APPEND(ols_check_input,__PREC) (y, x, beta, add_const, trans_x, &
        rcond, res, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: y
    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:,:) :: beta
    logical, intent(in) :: add_const
    logical, intent(in) :: trans_x
    real (PREC), intent(in) :: rcond
    type (__APPEND(lm_data,__PREC)), intent(in), dimension(:), optional :: res
    type (status_t), intent(out) :: status

    integer :: nvars, nobs, nconst, nlhs, ncoefs

    status = NF_STATUS_INVALID_ARG

    call ols_get_dims (y, x, add_const, trans_x, nobs, nvars, nlhs, ncoefs, nconst)

    if (nobs < ncoefs) return
    if (size(y, 1) /= nobs) return
    if (size(beta, 2) /= nlhs) return
    if (size(beta, 1) /= ncoefs) return
    if (rcond < 0) return

    if (present(res)) then
        if (size(res) /= nlhs) return
    end if

    status = NF_STATUS_OK

end subroutine



pure subroutine __APPEND(ols_get_dims,__PREC) (y, x, add_const, trans_x, &
        nobs, nvars, nlhs, ncoefs, nconst)
    !*  OLS_GET_DIMS returns the number of various dimensions of the least-squares
    !   problem, depending on input data arrays and whether these should be
    !   transposed, an intercept added, etc.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    logical, intent(in) :: add_const
    logical, intent(in) :: trans_x
    integer, intent(out) :: nobs, nvars, nlhs, ncoefs, nconst

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
    ncoefs = nvars + nconst

end subroutine



subroutine __APPEND(ols_2d,__PREC) (y, x, beta, add_const, trans_x, rcond, &
        rank, res, status)
    !*  OLS_2D computes the ordinary least-squares problem for given independent
    !   data X and (potentially multiple) dependent variables Y.
    !
    !   Note that a regression constant is added by default.
    !
    !   The LS problem is solved using SVD as implemented in LAPACK's GELSD
    !   routine. The optional arguments RCOND and RANK are passed directly
    !   to/from GELSD.
    !
    !   Note that the routine creates copies for input arrays X and Y as these
    !   will be overwritten by GELSD.

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: y
        !*  Array of LHS variables (separate regression is performed for each
        !   LHS variables using the same set of RHS variables)
    real (PREC), intent(in), dimension(:,:) :: x
        !*  Array of RHS variables
    real (PREC), intent(out), dimension(:,:) :: beta
        !*  Array of estimated coefficients
    logical, intent(in), optional :: add_const
        !*  If present and .TRUE., add constant to RHS variables (default: .TRUE.)
    logical, intent(in), optional :: trans_x
        !*  If present and .TRUE., transpose array X of RHS variables.
    real (PREC), intent(in), optional :: rcond
        !*  Optional argument RCOND passed to LAPACK's GELSD that allows to
        !   control the effective rank of the regressor matrix.
    integer, intent(out), optional :: rank
        !*  Contains effective rank of regressor matrix
    type (__APPEND(lm_data,__PREC)), dimension(:), optional :: res
        !*  Optional result objects for linear models. Note that a separate
        !   object is returned for each LHS variable.
    type (status_t), intent(out), optional :: status
        !*  Optional exit code

    logical :: ladd_const, ltrans_x
    real (PREC) :: lrcond, var_rhs
    integer :: nobs, nvars, ncoefs, nconst, nlhs, i
    type (status_t) :: lstatus

    ! GELSD arguments
    real (PREC), dimension(:), allocatable :: work, sval
    integer, dimension(:), allocatable :: iwork
    real (PREC), dimension(:,:), allocatable :: lx, lhs
    integer :: lrank
    integer :: lwork, info, m, n, nrhs, lda, ldb, mn, liwork

    lstatus = NF_STATUS_OK

    ladd_const = .true.
    if (present(add_const)) ladd_const = add_const

    ltrans_x = .false.
    if (present(trans_x)) ltrans_x = trans_x

    ! Use default value from LAPACK95 GELS wrapper
    lrcond = 100 * epsilon(1.0_PREC)
    if (present(rcond)) lrcond = rcond

    call ols_check_input (y, x, beta, ladd_const, ltrans_x, lrcond, res, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    call ols_get_dims (y, x, ladd_const, ltrans_x, nobs, nvars, nlhs, ncoefs, nconst)

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
    allocate (lhs(nobs,nlhs), source=y)

    ! 1. Workspace query
    lwork = -1
    allocate (work(1))
    allocate (iwork(1))

    call LAPACK_GELSD (m, n, nrhs, lx, lda, lhs, ldb, sval, lrcond, lrank, work, &
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
    call LAPACK_GELSD (m, n, nrhs, lx, lda, lhs, ldb, sval, lrcond, lrank, work, &
        lwork, iwork, info)

    ! Check whether algorithm for computing SVD failed to converge
    if (info > 0) then
        lstatus = NF_STATUS_NOT_CONVERGED
        goto 100
    end if

    ! Copy coefficients
    forall (i=1:nlhs) beta(1:ncoefs,i) = lhs(1:ncoefs,i)
    ! Copy over optional output arguments
    if (present(rank)) rank = lrank

    if (present(res)) then
        ! Update LM_DATA result objects for OLS model

        ! Fraction of RHS variance used, analogous to PCR regression
        ! (only applicable if RHS matrix does not have full rank)
        var_rhs = sum(sval(1:lrank) ** 2.0_PREC) / sum(sval**2.0_PREC)

        do i = 1, nlhs
            call lm_data_update (res(i), model=NF_STATS_LM_OLS, &
                coefs=beta(:,i), nobs=nobs, nvars=nvars, var_expl=var_rhs, &
                rank_rhs=lrank, trans_rhs=ltrans_x, add_const=ladd_const)
        end do
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine __APPEND(ols_1d,__PREC) (y, x, beta, add_const, trans_x, rcond, &
        rank, res, status)
    !*  OLS_1D provides a convenient wrapper for OLS_2D for one-dimensional
    !   input data (ie for regressions with a single dependent variable).
    !
    !   See the documentation for OLS_2D for details.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), target :: y
    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(out), dimension(:) :: beta
    logical, intent(in), optional :: add_const
    logical, intent(in), optional :: trans_x
    real (PREC), intent(in), optional :: rcond
    integer, intent(out), optional :: rank
    type (__APPEND(lm_data,__PREC)), intent(out), optional :: res
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer :: ptr_y
    real (PREC), dimension(:,:), allocatable :: beta2d
    integer :: nobs, ncoefs
    type (status_t) :: lstatus
    type (__APPEND(lm_data,__PREC)), dimension(1) :: res1d

    nobs = size(y)
    ncoefs = size(beta)

    ptr_y(1:nobs,1:1) => y
    allocate (beta2d(ncoefs,1), source=0.0_PREC)

    if (present(res)) then
        call ols (ptr_y, x, beta2d, add_const, trans_x, rcond, rank, &
            res=res1d, status=lstatus)
    else
        call ols (ptr_y, x, beta2d, add_const, trans_x, rcond, rank, &
            status=lstatus)
    end if

    ! Leave output array unmodified if OLS could not be performed correctly
    ! to be consistent with behavior of OLS_2D.
    if (lstatus == NF_STATUS_OK) then
        beta(:) = beta2d(:,1)
        if (present(res)) res = res1d(1)
    end if

    if (present(status)) status = lstatus

end subroutine



!-------------------------------------------------------------------------------
! PCA (Principal component analysis)

pure subroutine __APPEND(pca_check_input,__PREC) (x, scores, ncomp, trans_x, s, &
        loadings, mean_x, std_x, propvar, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:,:) :: scores
    integer, intent(in) :: ncomp
    logical, intent(in) :: trans_x
    real (PREC), intent(in), dimension(:), optional :: s
    real (PREC), intent(in), dimension(:,:), optional :: loadings
    real (PREC), intent(in), dimension(:), optional :: mean_x, std_x
    real (PREC), intent(in), dimension(:), optional :: propvar
    type (status_t), intent(out) :: status

    integer :: nobs, nvars

    status = NF_STATUS_INVALID_ARG

    call pca_get_dims (x, trans_x, nobs, nvars)

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


pure subroutine __APPEND(pca_get_dims,__PREC) (x, trans_x, nobs, nvars)
    !*  PCA_GET_DIMS returns the number of various dimensions of the PCA
    !   problem derived from input data and options (e.g. whether data should
    !   be transposed).

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: x
    logical, intent(in) :: trans_x
    integer, intent(out) :: nobs, nvars

    if (trans_x) then
        nobs = size(x, 2)
        nvars = size(x, 1)
    else
        nobs = size(x, 1)
        nvars = size(x, 2)
    end if

end subroutine


subroutine __APPEND(pca,__PREC) (x, scores, ncomp, center, scale, trans_x, &
        sval, loadings, mean_x, std_x, propvar, status)

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(out), dimension(:,:), contiguous :: scores
    integer, intent(in) :: ncomp
    logical, intent(in), optional :: center
    logical, intent(in), optional :: scale
    logical, intent(in), optional :: trans_x
    real (PREC), intent(out), dimension(:), optional :: sval
    real (PREC), intent(out), dimension(:,:), optional :: loadings
    real (PREC), intent(out), dimension(:), optional :: mean_x
    real (PREC), intent(out), dimension(:), optional :: std_x
    real (PREC), intent(out), dimension(:), optional :: propvar
    type (status_t), intent(out), optional :: status

    logical :: lscale, lcenter, ltrans_x
    type (status_t) :: lstatus
    real (PREC), dimension(:,:), allocatable :: x_n, x_tmp
    real (PREC), dimension(:), allocatable :: work
    integer :: nvars, nobs, i
    real (PREC), dimension(:), allocatable :: lmean_x, lstd_x

    ! Argument for GESDD
    real (PREC), dimension(:), allocatable :: ls
    real (PREC), dimension(:,:), allocatable :: vt
    real (PREC), dimension(1) :: qwork
    integer, dimension(:), allocatable :: iwork
    character (1), parameter :: jobz = 'O'
    real (PREC), dimension(0,0) :: u
    integer :: lda, ldu, ldvt, lwork, info, m, n, mn

    ! Arguments for GEMM
    real (PREC), dimension(:,:), allocatable :: vk
    character (1) :: transa, transb
    real (PREC), parameter :: alpha = 1.0_PREC, beta = 0.0_PREC
    integer :: ldb, ldc, k

    lstatus = NF_STATUS_OK

    lcenter = .true.
    lscale = .true.
    ltrans_x = .false.
    if (present(center)) lcenter = center
    if (present(scale)) lscale = scale
    if (present(trans_x)) ltrans_x = trans_x

    call pca_check_input (x, scores, ncomp, ltrans_x, sval, loadings, &
        mean_x, std_x, propvar, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    call pca_get_dims (x, ltrans_x, nobs, nvars)

    ! 1. Normalize data in place; hence we first need to create a copy,
    ! accounting for the possibility that x needs to be transformed.
    ! All code below assumed that x_n has shape (noabs, nvars).
    if (ltrans_x) then
        allocate (x_n(nobs, nvars))
        forall (i=1:nvars) x_n(:,i) = x(i,1:nobs)
    else
        allocate (x_n(nobs, nvars), source=x)
    end if

    allocate (lmean_x(nvars), lstd_x(nvars))

    if (lcenter .or. lscale) then
        ! Compute mean and std. dev. for each variable
        call std (x_n, s=lstd_x, m=lmean_x, dof=0, dim=1, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! de-mean variables
        if (lcenter) then
            do i = 1, nvars
                x_n(:,i) = x_n(:,i) - lmean_x(i)
            end do
        end if

        ! rescale variables to have std. dev. of 1
        if (lscale) then
            do i = 1, nvars
                ! Skip constant variables to avoid generating NaNs
                if (lstd_x(i) > 0.0_PREC) then
                    x_n(:,i) = x_n(:,i) / lstd_x(i)
                end if
            end do
        end if
    end if

    ! Notes on calling GESDD: we are using JOBZ='O' so that only the first
    ! NVARS columns of the matrix U in the SVD X=USV' will be computed and
    ! stored in the array X_TMP.
    ! The input checks ensure that NOBS >= NVARS, so the matrix V' will always
    ! be written into the array VT and never into X_TMP.
    ! For JOBZ='O' and m >= n, the array U will never be referenced, so we
    ! can just pass a size-0 dummy array.
    m = nobs
    n = nvars
    lda = nobs
    ldvt = n
    ldu = m
    mn = max(1, min(m, n))

    allocate (iwork(8*mn))
    allocate (vt(n, n))
    allocate (ls(mn))
    ! Contents of x_n will be overwritten by GESDD, make copy
    allocate (x_tmp(nobs,nvars), source=x_n)

    ! workspace query
    lwork = -1
    call LAPACK_GESDD (jobz, m, n, x_tmp, lda, ls, u, ldu, vt, ldvt, qwork, &
        lwork, iwork, info)

    ! Recover minimal work space size
    if (info /= 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if
    lwork = int(qwork(1))

    ! perform actual SVD
    allocate (work(lwork))

    call LAPACK_GESDD (jobz, m, n, x_tmp, lda, ls, u, ldu, vt, ldvt, work, &
        lwork, iwork, info)
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
    call BLAS_GEMM (transa, transb, m, n, k, alpha, x_n, lda, vk, ldb, beta, scores, ldc)

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

    integer :: nobs, nvars, ncomp, nlhs, nconst, ncoefs

    status = NF_STATUS_INVALID_ARG

    call pcr_get_dims (lhs, scores, coefs, add_const, &
        nobs, nvars, ncomp, nlhs, ncoefs, nconst)

    if (nobs < 1 .or. nobs < ncomp) return
    ! Check shapes for all arrays, even though some of these were used to
    ! obtain the dimensions in the first place and are thus by definition true.
    if (size(lhs,1) /= nobs .or. size(lhs,2) /= nlhs) return
    if (size(scores,1) /= nobs .or. size(scores,2) /= ncomp) return
    if (size(sval) < ncomp) return
    if (size(loadings,1) /= nvars .or. size(loadings,2) /= ncomp) return
    if (size(coefs,1) /= ncoefs .or. size(coefs,2) /= nlhs) return
    if (present(mean_x)) then
        if (size(mean_x) < nvars) return
    end if
    if (present(std_x)) then
        if (size(std_x) < nvars) return
    end if

    ! Do not support regressing on zero PCs if constant adding constant
    ! was not requested by user.
    if (ncomp == 0.0 .and. .not. add_const) return

    status = NF_STATUS_OK

end subroutine

pure subroutine __APPEND(pcr_get_dims,__PREC) (lhs, scores, coefs, add_const, &
        nobs, nvars, ncomp, nlhs, ncoefs, nconst)
    !*  PCR_GET_DIMS returns the number of various dimensions of the PCR
    !   problem, depending on input data arrays and whether an intercept
    !   shoud be added, etc.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: lhs, scores, coefs
    logical, intent(in) :: add_const
    integer, intent(out) :: nobs, nvars, nlhs, ncoefs, nconst, ncomp

    nconst = 0
    if (add_const) nconst = 1
    ncoefs = size(coefs, 1)
    nvars = ncoefs - nconst
    ncomp = size(scores, 2)
    nobs = size(lhs, 1)
    nlhs = size(lhs, 2)

end subroutine



subroutine __APPEND(pcr_2d,__PREC) (lhs, scores, sval, loadings, coefs, mean_x, std_x, &
        center, add_const, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: lhs
    real (PREC), intent(in), dimension(:,:), contiguous :: scores
    real (PREC), intent(in), dimension(:) :: sval
    real (PREC), intent(in), dimension(:,:), contiguous :: loadings
    real (PREC), intent(out), dimension(:,:) :: coefs
    real (PREC), intent(in), dimension(:), contiguous, optional :: mean_x
    real (PREC), intent(in), dimension(:), contiguous, optional :: std_x
    logical, intent(in), optional :: add_const
    logical, intent(in), optional :: center
    type (status_t), intent(out), optional :: status

    logical :: ladd_const, lcenter
    integer :: nvars, nobs, ncomp, nlhs, ncoefs
    real (PREC), dimension(:), pointer :: ptr_mean_x, ptr_std_x

    type (status_t) :: lstatus
    integer :: i, nconst

    real (PREC), dimension(:,:), allocatable :: y_n, PCy, lcoefs
    real (PREC), dimension(:), allocatable :: work
    real (PREC), dimension(:), allocatable :: mean_y
    ! Arguments to GEMM
    character (1) :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real (PREC) :: alpha, beta
    ! Arguments to GEMV
    integer, parameter :: incx = 1, incy = 1

    nullify (ptr_mean_x, ptr_std_x)

    lstatus = NF_STATUS_OK

    ladd_const = .true.
    lcenter = .true.
    if (present(add_const)) ladd_const = add_const
    if (present(center)) lcenter = center

    call pcr_check_input (lhs, scores, sval, loadings, coefs, mean_x, std_x, &
        ladd_const, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    call pcr_get_dims (lhs, scores, coefs, ladd_const, &
        nobs, nvars, ncomp, nlhs, ncoefs, nconst)

    ! copy over dependent variable, will be overwritten by GELS
    allocate (y_n(nobs, nlhs), source=lhs)
    allocate (mean_y(nlhs), source=0.0_PREC)
    if (lcenter) then
        call normalize (y_n, m=mean_y, dim=1, center=.true., scale=.false.)
    end if

    ! If no components are present, we can still run a constant-only regression
    ! if requested by user. Otherwise exit doing nothing.
    if (ncomp == 0) then
        if (ladd_const) then
            coefs(1,1:nlhs) = mean_y
        end if
        goto 100
    end if

    ! sample mean and standard deviation of original x variables
    call assert_alloc_ptr (mean_x, nvars, ptr_mean_x)
    call assert_alloc_ptr (std_x, nvars, ptr_std_x)

    if (.not. present(mean_x)) ptr_mean_x = 0.0_PREC
    if (.not. present(std_x)) ptr_std_x = 1.0_PREC

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
    alpha = 1.0_PREC
    beta = 0.0_PREC
    allocate (PCy(ncomp, nlhs))
    call BLAS_GEMM (transa, transb, m, n, k, alpha, scores, lda, y_n, ldb, beta, PCy, ldc)

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
    alpha = 1.0_PREC
    beta = 0.0_PREC
    allocate (lcoefs(nvars, nlhs))
    call BLAS_GEMM (transa, transb, m, n, k, alpha, loadings, lda, PCy, ldb, beta, lcoefs, ldc)

    ! Scale back betas if explanatory variables were normalized
    do i = 1, nlhs
        ! Prevent NaN's due to constant RHS variables
        allocate (work(size(std_x)), source=std_x)
        where (work == 0.0)
            work = 1.0
        end where

        lcoefs(:, i) = lcoefs(:, i) / work

        deallocate (work)
    end do

    ! Add intercept if requested
    if (ladd_const) then
        m = nvars
        n = nlhs
        lda = nvars
        transa = 'T'
        alpha = -1.0_PREC
        beta = 1.0

        ! Compute intercept as a = mean(y) - mean(x)'b
        ! where mean(x)'b is stored in temporary array work.
        ! Manually create temporary array as coefs(1:1,:) is not contiguous
        ! in memory and ifort issues runtime warning in debug mode every time.
        allocate (work(nlhs), source=mean_y)

        call BLAS_GEMV (transa, m, n, alpha, lcoefs, lda, mean_x, incx, beta, &
            work, incy)

        coefs(1,1:nlhs) = work

        deallocate (work)
    end if

    ! copy over remaining coefficients
    forall (i=1:nlhs) coefs(1+nconst:ncoefs,i) = lcoefs(:,i)

100 continue
    call assert_dealloc_ptr (mean_x, ptr_mean_x)
    call assert_dealloc_ptr (std_x, ptr_std_x)
    if (present(status)) status = lstatus

end subroutine



subroutine __APPEND(pcr_1d,__PREC) (lhs, scores, sval, loadings, coefs, &
        mean_x, std_x, add_const, center, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), target :: lhs
    real (PREC), intent(in), dimension(:,:), contiguous :: scores
    real (PREC), intent(in), dimensioN(:) :: sval
    real (PREC), intent(in), dimension(:,:), contiguous :: loadings
    real (PREC), intent(out), dimension(:), target :: coefs
    real (PREC), intent(in), dimension(:), contiguous, optional :: mean_x
    real (PREC), intent(in), dimension(:), contiguous, optional :: std_x
    logical, intent(in), optional :: add_const
    logical, intent(in), optional :: center
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer :: ptr_lhs, ptr_coefs

    integer :: ncoefs, nobs

    nobs = size(lhs)
    ncoefs = size(coefs)

    ptr_lhs(1:nobs,1:1) => lhs
    ptr_coefs(1:ncoefs,1:1) => coefs

    call pcr (ptr_lhs, scores, sval, loadings, ptr_coefs, &
        mean_x, std_x, add_const, center, status)

end subroutine


pure subroutine __APPEND(pcr_pca_get_dims,__PREC) (lhs, rhs, add_const, &
        trans_rhs, nobs, nvars, nlhs, ncoefs, nconst)
    !*  PCR_GET_DIMS returns the number of various dimensions of the PCR
    !   problem, depending on input data arrays and whether an intercept
    !   shoud be added, etc.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: lhs, rhs
    logical, intent(in) :: add_const
    logical, intent(in) :: trans_rhs
    integer, intent(out) :: nobs, nvars, nlhs, ncoefs, nconst

    nconst = 0
    if (add_const) nconst = 1

    if (trans_rhs) then
        nobs = size(rhs, 2)
        nvars = size(rhs, 1)
    else
        nobs = size(rhs, 1)
        nvars = size(rhs, 2)
    end if

    nlhs = size(lhs, 2)
    ncoefs = nvars + nconst

end subroutine


pure subroutine __APPEND(pcr_pca_check_input,__PREC) (lhs, rhs, ncomp_min, &
        coefs, add_const, trans_rhs, var_min, res, status)

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: lhs
    real (PREC), intent(in), dimension(:,:) :: rhs
    integer, intent(in), optional :: ncomp_min
    real (PREC), intent(in), dimension(:,:) :: coefs
    logical, intent(in) :: add_const
    logical, intent(in) :: trans_rhs
    real (PREC), intent(in), optional :: var_min
    type (__APPEND(lm_data,__PREC)), intent(in), dimension(:), optional :: res
    type (status_t), intent(out) :: status

    integer :: nobs, nvars, nlhs, nconst, ncoefs

    status = NF_STATUS_INVALID_ARG

    call pcr_pca_get_dims (lhs, rhs, add_const, trans_rhs, &
        nobs, nvars, nlhs, ncoefs, nconst)

    if (present(ncomp_min)) then
        if (ncomp_min < 1) return
        if (ncomp_min > nvars + nconst) return
        if (nobs < ncomp_min) return
    end if

    if (size(lhs,1) /= nobs .or. size(lhs,2) /= nlhs) return
    ! Check that coefficient array can hold coefs for all components
    ! Allow for more columns to be present, but not for more rows.
    if (size(coefs, 2) < nlhs) return
    if (ncoefs /= (nvars + nconst)) return
    if (present(var_min)) then
        if (var_min < 0.0_PREC .or. var_min > 1.0_PREC) return
    end if

    ! Make sure there is a LM_DATA object for each LHS variable
    if (present(res)) then
        if (size(res) /= nlhs) return
    end if

    status = NF_STATUS_OK

end subroutine


subroutine __APPEND(pcr_pca_2d,__PREC) (lhs, rhs, coefs, ncomp_min, add_const, &
        center, trans_rhs, var_min, res, status)

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: lhs
    real (PREC), intent(in), dimension(:,:) :: rhs
    real (PREC), intent(out), dimension(:,:) :: coefs
    integer, intent(in), optional :: ncomp_min
        !*  If present, specifies the min. number of principal components to
        !   use (default: all). If VAR_SHARE is not given, NCOMP_MIN will
        !   also be the actual number of components used.
    logical, intent(in), optional :: add_const
        !*  If present and true, an intercept will be automatically added
        !   to the RHS variables (default: true).
    logical, intent(in), optional :: center
        !*  If present and true, variables will be centered before running PCR
        !   (default: true).
    logical, intent(in), optional :: trans_rhs
        !*  If present and true, array of RHS variables will be transposed
        !   prior to running PCR (default: false).
    real (PREC), intent(in), optional :: var_min
        !*  If present, specifies the min. share of variance in the RHS
        !   variables that the (automatically) selected number of principal
        !   components should explain.
    type (__APPEND(lm_data,__PREC)), dimension(:), optional :: res
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    real (PREC), dimension(:,:), allocatable :: scores, loadings
    real (PREC), dimension(:), allocatable :: sval, mean_x, std_x, propvar
    integer :: nvars, nconst, nobs, ncoefs, nlhs, ncomp
    logical :: ltrans_rhs, lcenter, ladd_const
    real (PREC) :: var_expl
    integer :: k

    lstatus = NF_STATUS_OK

    ltrans_rhs = .false.
    lcenter = .true.
    ladd_const = .true.
    if (present(add_const)) ladd_const = add_const
    if (present(center)) lcenter = center
    if (present(trans_rhs)) ltrans_rhs = trans_rhs

    call pcr_pca_check_input (lhs, rhs, ncomp_min, coefs, ladd_const, &
        ltrans_rhs, var_min, res, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    call pcr_pca_get_dims (lhs, rhs, ladd_const, ltrans_rhs, &
        nobs, nvars, nlhs, ncoefs, nconst)

    if (.not. present(var_min) .and. present(ncomp_min)) then
        ncomp = ncomp_min
    else
        ncomp = nvars
    end if

    ! Allocate working arrays
    allocate (scores(nobs, ncomp))
    allocate (loadings(nvars, ncomp))
    allocate (sval(ncomp))
    allocate (mean_x(nvars), std_x(nvars))
    allocate (propvar(nvars), source=0.0_PREC)

    ! Perform principal component analysis
    call pca (rhs, scores, ncomp, center=.true., scale=.true., trans_x=ltrans_rhs, &
        sval=sval, loadings=loadings, mean_x=mean_x, std_x=std_x, &
        propvar=propvar, status=lstatus)
    if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

    ! Select number of principal components to use, if applicable
    if (present(var_min)) then
        var_expl = 0.0
        do k = 1, nvars
            var_expl = var_expl + propvar(k)
            if (var_expl >= var_min) exit
        end do
        ncomp = min(k, nvars)
        if (present(ncomp_min)) then
            ncomp = max(ncomp_min, ncomp)
        end if
    end if

    ! Run principal component regression
    call pcr (lhs, scores(:,1:ncomp), sval(1:ncomp), loadings(:,1:ncomp), &
        coefs, mean_x, std_x, lcenter, ladd_const, lstatus)
    if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

    if (present(res)) then
        var_expl = sum(propvar(1:ncomp))
        ! Update one LM_DATA object for each LHS
        do k = 1, size(res)
            call lm_data_update (res(k), model=NF_STATS_LM_PCR, coefs=coefs(:,k), &
                nobs=nobs, nvars=nvars, ncomp=ncomp, var_expl=var_expl, &
                add_const=ladd_const, trans_rhs=ltrans_rhs)
        end do
    end if

100 continue

    if (allocated(scores)) deallocate (scores)
    if (allocated(loadings)) deallocate (loadings)
    if (allocated(sval)) deallocate (sval)
    if (allocated(mean_x)) deallocate (mean_x)
    if (allocated(std_x)) deallocate (std_x)
    if (allocated(propvar)) deallocate (propvar)

    if (present(status)) status = lstatus

end subroutine


subroutine __APPEND(pcr_pca_1d,__PREC) (lhs, rhs, coefs, ncomp_min, add_const, &
        center, trans_rhs, var_min, res, status)

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), contiguous, target :: lhs
    real (PREC), intent(in), dimension(:,:) :: rhs
    real (PREC), intent(out), dimension(:), optional :: coefs
    integer, intent(in), optional :: ncomp_min
    logical, intent(in), optional :: add_const
    logical, intent(in), optional :: center
    logical, intent(in), optional :: trans_rhs
    real (PREC), intent(in), optional :: var_min
    type (__APPEND(lm_data,__PREC)), intent(in out), optional :: res
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_lhs
    real (PREC), dimension(:,:), allocatable :: coefs2d
    type (__APPEND(lm_data,__PREC)), dimension(1) :: res1d

    integer :: nobs, ncoefs, nconst, nvars, nlhs
    logical :: ladd_const, ltrans_rhs

    ladd_const = .false.
    ltrans_rhs = .false.
    if (present(add_const)) ladd_const = add_const
    if (present(trans_rhs)) ltrans_rhs = trans_rhs

    nobs = size(lhs)
    ptr_lhs(1:nobs,1:1) => lhs

    call pcr_pca_get_dims (ptr_lhs, rhs, ladd_const, ltrans_rhs, &
        nobs, nvars, nlhs, ncoefs, nconst)

    allocate (coefs2d(ncoefs,1))

    if (present(res)) then
        call pcr (ptr_lhs, rhs, coefs2d, ncomp_min, add_const, center, &
            trans_rhs, var_min, res1d, status)
        res = res1d(1)
    else
        call pcr (ptr_lhs, rhs, coefs2d, ncomp_min, add_const, center, &
            trans_rhs, var_min, status=status)
    end if

    if (present(coefs)) coefs = coefs2d(:,1)

end subroutine


subroutine __APPEND(lm_post_estim,__PREC) (model, lhs, rhs, rsq, status)
    !*  LM_POST_ESTIM computes common post-estimation statistics such as
    !   R^2 for a given linear model. The model must have been estimated
    !   using before invoking this routine.
    integer, parameter :: PREC = __PREC
    type (__APPEND(lm_data,__PREC)), intent(in) :: model
        !*  Estimated model
    real (PREC), intent(in), dimension(:), contiguous :: lhs
        !*  Vector of observations of the dependent variable used to estimate
        !   model.
    real (PREC), intent(in), dimension(:,:), contiguous :: rhs
        !*  Array of explanatory (RHS) variables used to estimate model.
        !   Array is assumed to have shape [NOBS, NVARS], unless
        !   TRANS_RHS is present and has value .TRUE.
    real (PREC), intent(out), optional :: rsq
        !*  If present, contains the model's R^2 on exit.
    type (status_t), intent(out), optional :: status

    logical :: need_resid
    real (PREC), dimension(:), allocatable :: resid

    type (status_t) :: lstatus
    real (PREC) :: var_u, var_y
    real (PREC), dimension(:), allocatable :: coefs
    integer :: ncoefs

    ! Variables used for BLAS routines
    integer, parameter :: incx = 1, incy = 1
    character (1) :: trans
    real (PREC) :: alpha, beta
    integer :: m, n, lda

    lstatus = NF_STATUS_OK

    ! List of outputs that require predicted values to compute
    need_resid = present(rsq)

    ! Check whether model has been estimated
    if (.not. allocated(model%coefs)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ncoefs = size(model%coefs)

    if (need_resid) then
        ! Compute YHAT using GEMV BLAS routine
        allocate (resid(model%nobs), source=lhs)

        if (model%add_const) then
            ! Constant was added during estimation process but is not present
            ! in RHS array
            allocate (coefs(model%nvars), source=model%coefs(2:ncoefs))
            ! Adjust LHS variable to compensate for "lacking" intercept in YHAT
            ! Note: At this point RESID contains just Y
            ! We want resid = y - yhat = y - (alpha + x'b) = (y-alpha) - x'b
            resid(:) = resid - model%coefs(1)
        else
            allocate (coefs(model%nvars), source=model%coefs)
        end if

        if (model%trans_rhs) then
            trans = 'T'
            m = model%nvars
            n = model%nobs
        else
            trans = 'N'
            m = model%nobs
            n = model%nvars
        end if

        lda = m
        ! Set up ALPHA, BETA such that result of GEMV will be resid = -yhat + y
        alpha = -1.0_PREC
        beta = 1.0_PREC
        call BLAS_GEMV (trans, m, n, alpha, rhs, lda, coefs, incx, beta, resid, incy)
        deallocate (coefs)
    end if

    if (present(rsq)) then
        ! Compute variance of residuals
        var_u = BLAS_DOT (model%nobs, resid, incx, resid, incy) / model%nobs

        ! Compute variable of LHS variable Y
        call std (lhs, s=var_y, dof=0)
        var_y = var_y ** 2.0_PREC

        rsq = 0.0_PREC
        if (var_y > 0.0_PREC) then
            rsq = 1.0_PREC - var_u / var_y
        end if
    end if


100 continue

    if (allocated(resid)) deallocate (resid)
    if (present(status)) status  = lstatus

end subroutine
