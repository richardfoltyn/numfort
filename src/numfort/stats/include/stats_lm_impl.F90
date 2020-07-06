

subroutine lm_result_update (self, conf, model, coefs, nobs, nrhs, nlhs, &
        ncomp, var_rhs, rank_rhs)

    type (lm_result), intent(inout) :: self
    integer, intent(in), optional :: model
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:), optional :: coefs
    integer, intent(in), optional :: nobs
    integer, intent(in), optional :: nrhs
    integer, intent(in), optional :: nlhs
    integer, intent(in), optional :: ncomp
    real (PREC), intent(in), optional :: var_rhs
    integer, intent(in), optional :: rank_rhs

    if (present(coefs)) then
        call copy_alloc (coefs, self%coefs)
    end if

    if (present(model)) self%model = model
    if (present(conf)) self%conf = conf
    if (present(nobs)) self%nobs = nobs
    if (present(nrhs)) self%nrhs = nrhs
    if (present(nlhs)) self%nlhs = nlhs
    if (present(ncomp)) self%ncomp = ncomp
    if (present(var_rhs)) self%var_rhs = var_rhs
    if (present(rank_rhs)) self%rank_rhs = rank_rhs
end subroutine



subroutine lm_result_reset (self)
    type (lm_result), intent(out) :: self

    ! INTENT(OUT) makes sure that all allocatable attribute arrays are
    ! deallocated.

    ! Overwrite with default LM_CONFIG, whatever that may be
    type (lm_config) :: conf
    self%conf = conf

    self%intercept = 0.0

    self%model = 0

    self%nrhs = 0
    self%nlhs = 0
    self%nobs = 0
    self%ncomp = 0

    self%rank_rhs = 0
    self%var_rhs = 0.0

end subroutine



pure subroutine get_dims_2d (X, Y, coefs, intercept, add_intercept, trans_x, &
        trans_y, trans_coefs, nobs, nrhs, nlhs, ncoefs, nconst_add, nconst_coefs, &
        status)
    !*  GET_DIMS returns the dimensions corresponding to various arrays
    !   used in linear models.
    real (PREC), intent(in), dimension(:,:) :: X
    real (PREC), intent(in), dimension(:,:) :: Y
    real (PREC), intent(in), dimension(:,:), optional :: coefs
    real (PREC), intent(in), dimension(:), optional :: intercept
    logical, intent(in), optional :: add_intercept
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
    logical, intent(in), optional :: trans_coefs
    integer, intent(out), optional :: nobs, nrhs, nlhs
    integer, intent(out), optional :: ncoefs
        !*  Number of coefficients stored in COEFS array.
    integer, intent(out), optional :: nconst_add
        !*  Number of constants ADDED to the model (irrespective of what is
        !   in the user-provided regressor matrix). Either 0 or 1, the latter
        !   being the case when ADD_INTERCEPT=.TRUE.
    integer, intent(out), optional :: nconst_coefs
        !*  Number of constants added to COEFS arrays. This value satisfies
        !   0 <= Nconst_coefs <= Nconst <= 1.
        !   Non-zero value can only be returned if COEFS is present, otherwise
        !   it is impossible to determine from the inputs.
    type (status_t), intent(out), optional :: status
        !*  Exit code. Returns NF_STATUS_INVALID_ARG if any of the dimensions
        !   are non-conformable with the given array sizes.

    logical :: ltrans_x, ltrans_y, ltrans_coefs
    logical :: dim_ok
    integer :: lnrhs, lnlhs, lncoefs, lnobs, lnobs_lhs, dim
    integer :: lnconst_add, lnconst_coefs
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'GET_DIMS'
    type (lm_config) :: conf

    lstatus = NF_STATUS_OK

    ! Use CONF values as default values
    call set_optional_arg (trans_x, conf%trans_x, ltrans_x)
    call set_optional_arg (trans_y, conf%trans_y, ltrans_y)
    call set_optional_arg (trans_coefs, conf%trans_coefs, ltrans_coefs)

    lnobs = 0
    lnrhs = 0
    lnlhs = 0
    lncoefs = 0
    lnconst_add = 0
    lnconst_coefs = 0

    if (present(add_intercept)) then
        if (add_intercept) lnconst_add = 1
    end if

    if (.not. ltrans_x) then
        lnobs = size(x, 1)
        lnrhs = size(x, 2)
    else
        lnobs = size(x, 2)
        lnrhs = size(x, 1)
    end if

    if (ltrans_y) then
        lnlhs = size(Y, 1)
        lnobs_lhs = size(Y, 2)
    else
        lnlhs = size(Y, 2)
        lnobs_lhs = size(Y, 1)
    end if

    if (present(coefs)) then
        ! Dimensions of COEFS array
        dim = 1
        if (ltrans_coefs) dim = 2

        lncoefs = size(coefs,dim)
        lnconst_coefs = lncoefs - lnrhs

        if (lnconst_coefs < 0 .or. lnconst_coefs > 1) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if

        if (present(add_intercept)) then
            if (.not. add_intercept) then
                ! Intercept should not be added, so don't allow COEFS
                ! to have one additional row/column
                dim_ok = (lncoefs == lnrhs) .and. (size(coefs,3-dim) == lnlhs)
            else
                ! Intercept to be added to model, but user still has the option
                ! of having it stored in COEFS or not, so both are acceptable
                dim_ok = ((lncoefs == lnrhs) .or. (lncoefs == (lnrhs + 1))) &
                    .and. (size(coefs,3-dim) == lnlhs)
            end if
        else
            ! ADD_INTERCEPT not present, so we don't know whether COEFS
            ! has one row/column more than NRHS.
            ! This is relevant when called from the PREDICT routines.
            dim_ok = (lncoefs == (lnrhs + 1) .or. lncoefs == lnrhs) &
                .and. size(coefs,3-dim) == lnlhs
        end if

        call check_cond (dim_ok, NAME, 'COEFS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    call check_cond (lnobs == lnobs_lhs, NAME, 'Y: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(intercept)) then
        if (size(intercept) < lnlhs) goto 100
    end if

100 continue

    if (present(nobs)) nobs = lnobs
    if (present(nrhs)) nrhs = lnrhs
    if (present(nlhs)) nlhs = lnlhs
    if (present(ncoefs)) ncoefs = lncoefs
    if (present(nconst_add)) nconst_add = lnconst_add
    if (present(nconst_coefs)) nconst_coefs = lnconst_coefs

    if (present(status)) status = lstatus

end subroutine



subroutine ols_2d (conf, X, Y, coefs, intercept, rank, res, status)
    !*  OLS_2D computes the ordinary least-squares problem for given independent
    !   data X and (potentially multiple) dependent variables Y.
    !
    !   The LS problem is solved using SVD as implemented in LAPACK's GELSD
    !   routine. The optional arguments RCOND and RANK are passed directly
    !   to/from GELSD.
    !
    !   Intercepts
    !   ----------
    !   The routine supports the following scenarios:
    !       1.  CONF%ADD_INTERCEPT = .FALSE. and INTERCEPT not present.
    !           No intercept should be added or returned: the model is
    !           estimated with an intercept if and only if the user
    !           manually adds a constant to X.
    !       2.  CONF%ADD_INTERCEPT = .TRUE. and INTERCEPT not present:
    !           The routine adds an intercept to the model.
    !       3.  CONF%ADD_INTERCEPT = .FALSE. and INTERCEPT is present:
    !           The routine adds an intercept to the model and its estimates
    !           are stored in INTERCEPT.
    !       4.  ADD_INTERCEPT = .TRUE. and INTERCEPT is present:
    !           The routine adds an intercept to the model and its estimates
    !           are stored INTERCEPT.
    !   Consequently, in cases (2) - (4), the regressor matrix is augmented
    !   to be [1, X], and the RANK output argument will be the rank of
    !   this matrix.
    !
    !   Additionally, if the array COEFS has shape [NRHS+1,NLHS], the
    !   estimated intercepts will ADDITIONALLY be stored in the first
    !   row of COEFS.
    !
    !   Note that the routine creates copies for input arrays X and Y as these
    !   will be overwritten by GELSD.
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  Array of RHS variables
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
        !*  Array of LHS variables (separate regression is performed for each
        !   LHS variables using the same set of RHS variables)
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: coefs
        !*  Array of estimated coefficients. Argument is optional in case
        !   coefficients should only be stored in RES argument.
    real (PREC), intent(out), dimension(:), optional :: intercept
    integer, intent(out), optional :: rank
        !*  Contains effective rank of regressor matrix
    type (lm_result), intent(inout), optional :: res
        !*  Optional result object for linear models.
    type (status_t), intent(out), optional :: status
        !*  Optional exit code

    real (PREC) :: var_rhs
    integer :: nobs, Nrhs, Ncoefs, nconst_add, nlhs, i
    type (status_t) :: lstatus
    type (lm_config) :: lconf
    character (*), parameter :: NAME = 'OLS'

    ! GELSD arguments
    real (PREC), dimension(:), allocatable :: work, sval
    integer, dimension(:), allocatable :: iwork
    real (PREC), dimension(:,:), allocatable :: Xp, Yp
    integer :: lrank
    integer :: lwork, info, m, n, lda, ldb, mn, liwork

    lstatus = NF_STATUS_OK

    ! --- Input processsing ---

    if (present(conf)) lconf = conf

    call check_cond (lconf%rcond >= 0.0_PREC, NAME, 'CONF: Invalid value for RCOND', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call get_dims (X, Y, coefs, intercept, lconf%add_intercept, lconf%trans_x, &
        lconf%trans_y, trans_coefs=lconf%trans_coefs, nobs=nobs, nrhs=Nrhs, &
        nlhs=Nlhs, ncoefs=Ncoefs, nconst_add=Nconst_add, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (Nobs >= (Nrhs + Nconst_add), NAME, 'Too few observations', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(intercept)) then
        call check_cond (size(intercept) >= nlhs, NAME, &
            'INTERCEPT: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

   ! --- Pre-process data ---

    ! Allocate (possibly transposed) array to store X variables; will be
    ! overwritten by GELSD, so we have to allocate in any case.
    ! Add constant as requested.
    allocate (Xp(nobs,Nrhs+Nconst_add))

    call transform_regressors (X, Xp, trans=lconf%trans_x, drop_na=.false., &
        add_intercept=lconf%add_intercept, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (lconf%trans_y) then
        allocate (Yp(nobs,Nlhs))
        call transform_regressors (Y, Yp, trans=lconf%trans_y, drop_na=.false., &
            status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        ! Just create a copy since Y would otherwise be overwritten in place
        allocate (Yp(nobs,nlhs), source=Y)
    end if

    ! --- GELSD ---

    m = nobs
    n = Nrhs + Nconst_add
    mn = max(1, min(m, n))
    lda = nobs
    ldb = nobs
    ! GELSD solves Ax=b, we solve Y = Xb, so our nlhs = nrhs in GELSD notation
    allocate (sval(mn))

    ! 1. Workspace query
    lwork = -1
    allocate (work(1))
    allocate (iwork(1))

    call LAPACK_GELSD (m, n, nlhs, Xp, lda, Yp, ldb, sval, lconf%rcond, lrank, &
        work, lwork, iwork, info)

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
    call LAPACK_GELSD (m, n, nlhs, Xp, lda, Yp, ldb, sval, lconf%rcond, lrank, &
        work, lwork, iwork, info)

    ! Check whether algorithm for computing SVD failed to converge
    if (info > 0) then
        lstatus = NF_STATUS_NOT_CONVERGED
        goto 100
    else if (info < 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! --- Store results ---

    ! Copy coefficients
    if (present(coefs)) then
        if (lconf%trans_coefs) then
            if (Ncoefs == Nrhs) then
                ! No intercept should be added to COEFS
                forall (i=1:Nrhs) coefs(1:Nlhs,i) = Yp(Nconst_add+i,1:Nlhs)
            else
                ! Include added intercept in COEFS
                forall (i=1:Ncoefs) coefs(1:Nlhs,i) = Yp(i,1:Nlhs)
            end if
        else
            if (Ncoefs == Nrhs) then
                ! No intercept should be added to COEFS
                forall (i=1:nlhs) coefs(1:Nrhs,i) = Yp(1+Nconst_add:Nrhs+Nconst_add,i)
            else
                ! Include added intercept in COEFS
                forall (i=1:nlhs) coefs(1:Ncoefs,i) = Yp(1:Ncoefs,i)
            end if
        end if
    end if

    if (present(intercept)) then
        ! Costant was added to X, intercept is in the first row of Yp
        intercept = 0.0
        if (lconf%add_intercept) then
            intercept(1:Nlhs) = Yp(1,:)
        end if
    end if

    ! Copy over optional output arguments
    if (present(rank)) rank = lrank

    if (present(res)) then
        ! Update LM_RESULT object for OLS model

        ! Fraction of RHS variance used, analogous to PCR regression
        ! (only applicable if RHS matrix does not have full rank)
        var_rhs = sum(sval(1:lrank) ** 2.0_PREC) / sum(sval**2.0_PREC)

        call lm_result_update (res, conf=lconf, model=NF_STATS_LM_OLS, &
            nobs=nobs, nrhs=Nrhs, nlhs=Nlhs, var_rhs=var_rhs, rank_rhs=lrank)

        ! Never store intercept in first row in RES%COEFS_MULTI array
        call copy_alloc (Yp(1+Nconst_add:Nrhs+Nconst_add,1:Nlhs), res%coefs_multi)

        if (lconf%add_intercept) then
            call copy_alloc (Yp(1,1:Nlhs), res%intercept_multi)
        else
            ! We cannot leave the array unallocated or initialized, make
            ! sure that we set intercepts to zero.
            call cond_alloc (res%intercept_multi, Nlhs)
            res%intercept_multi(:) = 0.0
        end if
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine ols_1d (conf, X, y, coefs, intercept, rank, res, status)
    !*  OLS_1D provides a convenient wrapper for OLS_2D for one-dimensional
    !   input data (ie for regressions with a single dependent variable).
    !
    !   See the documentation for OLS_2D for details.
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: coefs
    real (PREC), intent(out), optional :: intercept
    integer, intent(out), optional :: rank
    type (lm_result), intent(inout), optional :: res
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), contiguous, pointer :: ptr_y
    real (PREC), dimension(:,:), contiguous, pointer :: ptr_coefs
    real (PREC), dimension(:), allocatable :: lintercept
    type (status_t) :: lstatus
    type (lm_config) :: lconf

    nullify (ptr_y, ptr_coefs)

    ptr_y(1:size(y),1:1) => y

    if (present(conf)) then
        ! Enforce no transposing of data in 1d arrays
        lconf = conf
        lconf%trans_y = .false.
        lconf%trans_coefs = .false.
    end if

    if (present(coefs)) then
        ptr_coefs(1:size(coefs),1:1) => coefs
    end if

    if (present(intercept)) then
        allocate (lintercept(1))
    end if

    call ols (conf, X, ptr_y, ptr_coefs, lintercept, rank, res, lstatus)

    if (lstatus == NF_STATUS_OK) then
        if (present(res)) then
            if (allocated(res%coefs_multi)) then
                call copy_alloc (res%coefs_multi(:,1), res%coefs)
            end if

            res%intercept = 0.0
            if (allocated(res%intercept_multi)) then
                res%intercept = res%intercept_multi(1)
            end if
        end if
    end if

    if (present(intercept)) then
        intercept = lintercept(1)
    end if

    if (present(status)) status = lstatus

end subroutine



subroutine pca (X, trans_x, center, scale, scores, sval, loadings, shift_x, &
        scale_x, propvar, status)

    real (PREC), intent(in), dimension(:,:), contiguous :: X
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: center
    logical, intent(in), optional :: scale
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: scores
    real (PREC), intent(out), dimension(:), contiguous, optional :: sval
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: loadings
    real (PREC), intent(out), dimension(:), contiguous, optional :: shift_x
    real (PREC), intent(out), dimension(:), contiguous, optional :: scale_x
    real (PREC), intent(out), dimension(:), contiguous, optional :: propvar
    type (status_t), intent(out), optional :: status

    logical :: lscale, lcenter, ltrans_x
    type (status_t) :: lstatus
    real (PREC), dimension(:,:), allocatable :: Xp
    real (PREC), dimension(:), allocatable :: work
    integer :: Nvars, Nobs, i, Ncomp
    character (*), parameter :: NAME = 'PCA'

    ! Argument for GESDD
    real (PREC), dimension(:), allocatable :: ls
    real (PREC), dimension(:,:), allocatable :: vt
    real (PREC), dimension(1) :: qwork
    integer, dimension(:), allocatable :: iwork
    character (1), parameter :: jobz = 'O'
    real (PREC), dimension(0,0) :: u
    integer :: lda, ldu, ldvt, lwork, info, m, n, mn

    lstatus = NF_STATUS_OK

    ! --- Input processing ---
    call set_optional_arg (center, .true., lcenter)
    call set_optional_arg (scale, .true., lscale)
    call set_optional_arg (trans_x, .false., ltrans_x)

    call get_regressor_dims (X, ltrans_x, Nvars=Nvars, Nobs=Nobs)

    ! --- Quick termination ---

    call check_cond (Nobs >= Nvars, NAME, 'X: too few observations', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Check remaining inputs ---

    if (present(scores)) then
        ! Do not check second dimension, we will simply fill in as many scores
        ! as there are columns
        call check_cond (size(scores, 1) == Nobs, NAME, &
            'SCORES: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(loadings)) then
        ! Don't check second dimension, we will simply fill in as many
        ! principal components as there are columns
        call check_cond (size(loadings,1) == Nvars, NAME, &
            'LOADINGS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Pre-process data ---

    ! 1. Standardize data in place; hence we first need to create a copy,
    ! accounting for the possibility that X needs to be transformed.
    ! GESDD overwrites its inputs, so we need a copy of X in any case!

    allocate (Xp(Nobs,Nvars))

    call transform_regressors (X, Xp, trans=ltrans_x, center=lcenter, &
        scale=lscale, shift_x=shift_x, scale_x=scale_x, drop_na=.false., &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- GESDD call ---


    ! Notes on calling GESDD: we are using JOBZ='O' so that only the first
    ! NVARS columns of the matrix U in the SVD X=USV' will be computed and
    ! stored in the array X_TMP.
    ! The input checks ensure that NOBS >= NVARS, so the matrix V' will always
    ! be written into the array VT and never into X_TMP.
    ! For JOBZ='O' and m >= n, the array U will never be referenced, so we
    ! can just pass a size-0 dummy array.
    m = Nobs
    n = Nvars
    lda = Nobs
    ldvt = n
    ldu = m
    mn = max(1, min(m, n))

    ! These arrays need to be allocated irrespective of whether GESDD is
    ! called or not
    allocate (vt(n, n))
    allocate (ls(mn))

    ! Skip the remaining SVD decomposition of no variables are present
    if (Nvars == 0) goto 50

    allocate (iwork(8*mn))

    ! workspace query
    lwork = -1
    call LAPACK_GESDD (jobz, m, n, Xp, lda, ls, u, ldu, vt, ldvt, qwork, &
        lwork, iwork, info)

    ! Recover minimal work space size
    if (info /= 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if
    lwork = int(qwork(1))

    ! perform actual SVD
    allocate (work(lwork))

    call LAPACK_GESDD (jobz, m, n, Xp, lda, ls, u, ldu, vt, ldvt, work, &
        lwork, iwork, info)

    if (info /= 0) then
        if (info < 0) then
            ! This should not happen since we checked arguments before
            lstatus = NF_STATUS_INVALID_ARG
        else
            lstatus = NF_STATUS_NOT_CONVERGED
        end if
        goto 100
    end if

    ! --- Store results ---

    ! Jump to this point if SVD decomposition is skipped
50  continue

    if (present(scores)) then
        ! Compute scores as S(:,i) = U(:,i) * sval(i) = X * V(:,i)
        ncomp = min(Nvars, size(scores, 2))
        do i = 1, ncomp
            scores(1:Nobs,i) = Xp(1:Nobs,i) * ls(i)
        end do

        ! Set any excess columns to zero
        scores(:,ncomp+1:) = 0.0
    end if

    if (present(sval)) then
        sval = 0.0
        ncomp = min(Nvars, size(sval))
        sval(1:ncomp) = ls(1:ncomp)
    end if

    if (present(loadings)) then
        ncomp = min(Nvars, size(loadings, 2))
        ! V^T returned by GESDD has shape [NCOMP,NVARS] where in that routine
        ! NCOMP = NVARS. If we want to restrict to NCOMP < NVARS, we have
        ! to use the first 1:NCOMP rows.
        do i = 1, ncomp
            loadings(1:Nvars,i) = vt(i,1:Nvars)
        end do
        ! Set any excess columns to zero
        loadings(:,ncomp+1:) = 0.0
    end if

    ! compute vector of variance share captured by each PC.
    ! Variance share is given by relative size of each eigenvalue (ie. squared
    ! singular value) relative to sum of all EVs.
    if (present(propvar)) then
        propvar(:) = 0.0
        if (allocated(ls)) then
            ncomp = min(Nvars, size(propvar), size(ls))
            propvar(1:ncomp) = ls(1:ncomp) ** 2 / sum(ls(1:Nvars) ** 2)
        end if
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine pcr_get_dims (Y, coefs, scores, sval, loadings, mean_x, &
        scale_x, add_intercept, trans_y, trans_coefs, nobs, nrhs, &
        ncomp, nlhs, ncoefs, nconst_add, nconst_coefs, status)
    !*  PCR_GET_DIMS returns the number of various dimensions of the PCR
    !   problem, depending on input data arrays and whether an intercept
    !   shoud be added, etc.
    real (PREC), intent(in), dimension(:,:) :: Y
    real (PREC), intent(in), dimension(:,:) :: coefs
    real (PREC), intent(in), dimension(:,:) :: scores
    real (PREC), intent(in), dimension(:) :: sval
    real (PREC), intent(in), dimension(:,:) :: loadings
    real (PREC), intent(in), dimension(:), optional :: mean_x, scale_x
    logical, intent(in), optional :: add_intercept
    logical, intent(in) :: trans_y
    logical, intent(in) :: trans_coefs
    integer, intent(out), optional :: nobs, nrhs, nlhs, ncoefs, nconst_add, ncomp
    integer, intent(out), optional :: nconst_coefs
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: lnobs, lnrhs, lnlhs, lncoefs, lnconst_add, lncomp, dim
    integer :: lnconst_coefs
    logical :: dim_ok
    character (*), parameter :: NAME = 'PCR_GET_DIMS'

    lstatus = NF_STATUS_INVALID_ARG

    lnobs = 0
    lnrhs = 0
    lnlhs = 0
    lncoefs = 0
    lnconst_add = 0
    lnconst_coefs = 0
    lncomp = 0

    if (present(add_intercept)) then
        if (add_intercept) lnconst_add = 1
    end if

    dim = 1
    if (trans_y) dim = 2
    lnobs = size(Y, dim)
    lnlhs = size(Y, 3-dim)

    lnrhs = size(loadings, 1)
    lncomp = size(scores, 2)

    if (size(scores,1) /= lnobs) goto 100
    if (size(loadings,2) /= lncomp) goto 100

    ! --- Check COEFS arrays ---

    dim = 1
    if (trans_coefs) dim = 2
    lncoefs = size(coefs,dim)
    lnconst_coefs = lncoefs - lnrhs

    if (lnconst_coefs < 0 .or. lnconst_coefs > 1) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(add_intercept)) then
        if (.not. add_intercept) then
            ! Intercept should not be added, so don't allow COEFS
            ! to have one additional row/column
            dim_ok = (lncoefs == lnrhs) .and. (size(coefs,3-dim) == lnlhs)
        else
            ! Intercept to be added to model, but user still has the option
            ! of having it stored in COEFS or not, so both are acceptable
            dim_ok = ((lncoefs == lnrhs) .or. (lncoefs == (lnrhs + 1))) &
                .and. (size(coefs,3-dim) == lnlhs)
        end if
    else
        ! ADD_INTERCEPT not present, so we don't know whether COEFS
        ! has one row/column more than NRHS.
        ! This is relevant when called from the PREDICT routines.
        dim_ok = (lncoefs == (lnrhs + 1) .or. lncoefs == lnrhs) &
            .and. size(coefs,3-dim) == lnlhs
    end if

    call check_cond (dim_ok, NAME, 'COEFS: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(mean_x)) then
        if (size(mean_x) /= lnrhs) goto 100
    end if

    if (present(scale_x)) then
        if (size(scale_x) /= lnrhs) goto 100
    end if

    lstatus = NF_STATUS_OK

100 continue

    if (present(nobs)) nobs = lnobs
    if (present(nrhs)) nrhs = lnrhs
    if (present(nlhs)) nlhs = lnlhs
    if (present(ncoefs)) ncoefs = lncoefs
    if (present(nconst_add)) nconst_add = lnconst_add
    if (present(nconst_coefs)) nconst_coefs = lnconst_coefs
    if (present(ncomp)) ncomp = lncomp

    if (present(status)) status = lstatus

end subroutine




subroutine pcr_2d (Y, scores, sval, loadings, coefs, shift_x, scale_x, &
        add_intercept, trans_y, trans_coefs, center, intercept, status)
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
    real (PREC), intent(in), dimension(:,:), contiguous :: scores
    real (PREC), intent(in), dimension(:), contiguous :: sval
    real (PREC), intent(in), dimension(:,:), contiguous :: loadings
    real (PREC), intent(out), dimension(:,:), contiguous :: coefs
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: shift_x
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: scale_x
    logical, intent(in), optional :: add_intercept
    logical, intent(in), optional :: trans_y
    logical, intent(in), optional :: trans_coefs
    logical, intent(in), optional :: center
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: intercept
    type (status_t), intent(out), optional :: status

    logical :: lcenter, ltrans_y, ltrans_coefs
    integer :: Nrhs, Nobs, ncomp, Nlhs, Ncoefs, Nconst_coefs, dim
    real (PREC), dimension(:), pointer, contiguous :: ptr_intercept
    integer :: i
    type (status_t) :: lstatus

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Yp
    real (PREC), dimension(:,:), allocatable :: PCy, lcoefs
    real (PREC), dimension(:), allocatable :: shift_y
    ! Arguments to GEMM
    character (1) :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real (PREC) :: alpha, beta
    ! Arguments to GEMV

    nullify (ptr_intercept, ptr_Yp)

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (center, .true., lcenter)
    call set_optional_arg (trans_y, .false., ltrans_y)
    call set_optional_arg (trans_coefs, .false., ltrans_coefs)

    call pcr_get_dims (Y, coefs, scores, sval, loadings, shift_x, scale_x, &
        add_intercept, ltrans_y, ltrans_coefs, nobs=Nobs, nrhs=Nrhs, &
        ncomp=Ncomp, nlhs=Nlhs, ncoefs=Ncoefs, nconst_coefs=Nconst_coefs, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Pre-process data ---

    if (ltrans_y .or. lcenter) then
        allocate (ptr_Yp(Nobs,Nlhs))
        allocate (shift_y(Nlhs))
        call transform_regressors (Y, ptr_Yp, center=lcenter, scale=.false., &
            trans=ltrans_y, shift_x=shift_y, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        ptr_Yp => Y
    end if

    ! If no components are present, we can still run a constant-only regression
    ! if requested by user. Otherwise exit doing nothing.
    if (ncomp == 0) then
        allocate (lcoefs(Nrhs,Nlhs), source=0.0_PREC)
        goto 50
    end if

    ! --- Compute slope coefficients ---

    ! Compute regression coefficients of regression Y on PC
    ! 1. Compute PC'y
    m = ncomp
    n = Nlhs
    k = Nobs
    lda = Nobs
    ldb = Nobs
    ldc = ncomp
    transa = 'T'
    transb = 'N'
    alpha = 1.0_PREC
    beta = 0.0_PREC
    allocate (PCy(ncomp, Nlhs))
    call BLAS_GEMM (transa, transb, m, n, k, alpha, scores, lda, ptr_Yp, ldb, beta, PCy, ldc)

    ! 2. Regression coefficients are given by (diag(sval**2))^{-1} PC'y,
    ! can be computed directly without solving equation system.
    do i = 1, Nlhs
        PCy(1:ncomp, i) = PCy(1:ncomp, i) * (sval(1:ncomp) ** (-2.0_PREC))
    end do

    ! PCR coefficients derived from regressing on k PC are given by
    ! beta_pcr = V_k * beta_ols
    m = Nrhs
    n = Nlhs
    k = ncomp
    lda = Nrhs
    ldb = ncomp
    ldc = Nrhs
    transa = 'N'
    transb = 'N'
    alpha = 1.0_PREC
    beta = 0.0_PREC
    allocate (lcoefs(Nrhs, Nlhs))
    call BLAS_GEMM (transa, transb, m, n, k, alpha, loadings, lda, PCy, ldb, beta, lcoefs, ldc)

    ! Scale back betas if explanatory variables were normalized
    if (present(scale_x)) then
        do i = 1, Nlhs
            lcoefs(1:Nrhs, i) = lcoefs(1:Nrhs, i) / scale_x(1:Nrhs)
        end do
    end if

50 continue

    ! --- Recover intercept ---

    dim = 1
    if (ltrans_coefs) dim = 2

    if (present(intercept) .or. size(coefs,1) > Nrhs) then

        call assert_alloc_ptr (intercept, Nlhs, ptr_intercept)

        ptr_intercept(:) = 0.0

        if (allocated(shift_y)) then
            ! If LHS variables were centered in this routine, initialize
            ! intercept using their mean (ie the value by which they were
            ! shifted).
            ! SHIFT_Y will be zeros if Y is pre-centered by TRANS_Y=.TRUE.
            ptr_intercept(1:Nlhs) = shift_y
        end if

        if (present(shift_x) .and. ncomp > 0) then
            ! Compute intercept as a = mean(y) - mean(x)'b
            ! where mean(x)'b is stored in temporary array work.
            ! Manually create temporary array as coefs(1:1,:) is not contiguous
            ! in memory and ifort issues runtime warning in debug mode every time.

            call BLAS_GEMV ('T', Nrhs, Nlhs, -1.0_PREC, lcoefs, Nrhs, &
                shift_x, 1, 1.0_PREC, ptr_intercept, 1)
        end if

        if (Nconst_coefs == 1) then
            if (ltrans_coefs) then
                coefs(1:Nlhs,1) = ptr_intercept
            else
                coefs(1,1:Nlhs) = ptr_intercept
            end if
        end if
    end if

    ! --- Store output ---

    if (ltrans_coefs) then
        forall (i=1:Nrhs) coefs(1:Nlhs,Nconst_coefs+i) = lcoefs(i,1:Nlhs)
    else
        forall (i=1:Nlhs) coefs(1+Nconst_coefs:Ncoefs,i) = lcoefs(:,i)
    end if

100 continue

    call assert_dealloc_ptr (Y, ptr_Yp)
    call assert_dealloc_ptr (intercept, ptr_intercept)

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_1d (y, scores, sval, loadings, coefs, shift_x, scale_x, &
        add_intercept, center, intercept, status)

    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(in), dimension(:,:), contiguous :: scores
    real (PREC), intent(in), dimensioN(:), contiguous :: sval
    real (PREC), intent(in), dimension(:,:), contiguous :: loadings
    real (PREC), intent(out), dimension(:), contiguous, target :: coefs
    real (PREC), intent(in), dimension(:), contiguous, optional :: shift_x
    real (PREC), intent(in), dimension(:), contiguous, optional :: scale_x
    logical, intent(in), optional :: add_intercept
    logical, intent(in), optional :: center
    real (PREC), intent(out), optional :: intercept
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), contiguous, pointer :: ptr_y, ptr_coefs
    real (PREC), dimension(:), allocatable :: lintercept
    integer :: ncoefs, nobs
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    nobs = size(y)
    ncoefs = size(coefs)

    ptr_y(1:nobs,1:1) => y
    ptr_coefs(1:ncoefs,1:1) => coefs

    if (present(intercept)) then
        allocate (lintercept(1))
    end if

    call pcr (ptr_y, scores, sval, loadings, ptr_coefs, shift_x, scale_x, &
        add_intercept, center=center, intercept=lintercept, status=lstatus)

    if (lstatus == NF_STATUS_OK) then
        if (present(intercept)) then
            intercept = lintercept(1)
        end if
    end if

end subroutine




subroutine pcr_pca_2d (conf, X, Y, coefs, ncomp, intercept, res, status)
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
    real (PREC), intent(out), dimension(:,:), contiguous, optional, target :: coefs
    integer, intent(out), optional :: ncomp
        !   If present, contains the actual number of PCs used on exit.
    real (PREC), intent(out), dimension(:), optional, contiguous, target :: intercept
    type (lm_result), intent(inout), optional, target :: res
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (lm_config) :: lconf
    real (PREC), dimension(:,:), allocatable :: lcoefs
    real (PREC), dimension(:), allocatable :: lintercept
    real (PREC), dimension(:,:), allocatable :: scores, loadings
    real (PREC), dimension(:), allocatable :: sval, shift_x, scale_x, propvar
    integer :: Nrhs, Nconst_add, nobs, ncoefs, Nlhs, lncomp, Nconst_coefs
    logical :: trans_coefs
    real (PREC) :: var_rhs
    integer :: k, i
    character (*), parameter :: NAME = 'PCR_PCA'

    lstatus = NF_STATUS_OK

    if (present(conf)) lconf = conf

    ! --- Input checks ---

    call get_dims (X, Y, coefs, intercept, lconf%add_intercept, lconf%trans_x, &
        lconf%trans_y, lconf%trans_coefs, nobs=nobs, nrhs=Nrhs, nlhs=nlhs, &
        ncoefs=ncoefs, nconst_add=Nconst_add, nconst_coefs=Nconst_coefs, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (Nobs >= (Nrhs + Nconst_add), NAME, 'Too few observations', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lconf%var_rhs_min >= 0.0_PREC .and. &
        lconf%var_rhs_min <= 1.0_PREC, NAME, 'CONF: Invalid value for VAR_RHS_MIN', &
        lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lconf%ncomp_min >= 0 .and. lconf%ncomp_min <= Nrhs, &
        NAME, 'CONF: Invalid value for NCOMP_MIN', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (lconf%ncomp /= DEFAULT_PCR_NCOMP) then
        call check_cond (lconf%ncomp >= 0 .and. lconf%ncomp <= Nrhs, NAME, &
            'CONF: Invalid value for NCOMP', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Run PCA ---

    ! Allocate working arrays
    ! Preliminary number of components: use all available, restrict below.
    lncomp = Nrhs
    allocate (scores(nobs, lncomp))
    allocate (loadings(Nrhs, lncomp))
    allocate (sval(lncomp))
    allocate (shift_x(Nrhs), scale_x(Nrhs))
    allocate (propvar(Nrhs), source=0.0_PREC)

    ! Skip PCA if no RHS are present
    if (Nrhs == 0) goto 50

    ! Perform principal component analysis
    call pca (X, center=lconf%center, scale=lconf%scale, trans_x=lconf%trans_x, &
        scores=scores, sval=sval, loadings=loadings, shift_x=shift_x, &
        scale_x=scale_x, propvar=propvar, status=lstatus)
    if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

    ! Select number of principal components to use:
    ! 1. Select by min. variance, if applicable;
    var_rhs = 0.0
    do k = 1, Nrhs
        var_rhs = var_rhs + propvar(k)
        if ((var_rhs - lconf%var_rhs_min) >= 0.0_PREC) exit
    end do
    lncomp = min(k, Nrhs)

    ! 2. Select by NCOMP_MIN
    lncomp = max(lconf%ncomp_min, lncomp)

    ! 3. User-imposed fixed number of components
    if (lconf%ncomp >= 0 .and. lconf%ncomp <= Nrhs) then
        lncomp = lconf%ncomp
    end if

    ! Jump to this point of PCA is skipped because of no RHS variables.
50  continue

    allocate (lcoefs(Nrhs,Nlhs))
    allocate (lintercept(Nlhs))
    ! Save user-provided TRANS_COEFS
    trans_coefs = lconf%trans_coefs
    lconf%trans_coefs = .false.

    ! Run principal component regression
    call pcr (Y, scores(:,1:lncomp), sval(1:lncomp), loadings(:,1:lncomp), &
        lcoefs, shift_x, scale_x, lconf%add_intercept, lconf%trans_y,  &
        lconf%trans_coefs, lconf%center, lintercept, lstatus)

    ! Restore user-provided TRANS_COEFS
    lconf%trans_coefs = trans_coefs

    if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

    ! --- Store result object ---

    if (present(intercept)) then
        intercept = 0.0
        intercept(1:Nlhs) = lintercept(1:Nlhs)
    end if

    if (present(coefs)) then
        if (lconf%trans_coefs) then
            do i = 1, Nlhs
                coefs(i,1+Nconst_coefs:Ncoefs) = lcoefs(1:Nrhs,i)
            end do
            if (Nconst_coefs == 1) then
                coefs(1:Nlhs,1) = lintercept
            end if
        else
            do i = 1, Nlhs
                coefs(1+Nconst_coefs:Ncoefs,i) = lcoefs(1:Nrhs,i)
            end do
            if (Nconst_coefs == 1) then
                coefs(1,1:Nlhs) = lintercept
            end if
        end if
    end if

    if (present(res)) then
        var_rhs = sum(propvar(1:lncomp))
        call lm_result_update (res, lconf, model=NF_STATS_LM_PCR, &
            nobs=nobs, nrhs=Nrhs, nlhs=Nlhs, ncomp=lncomp, var_rhs=var_rhs)

        if (allocated(res%coefs_multi)) then
            call copy_alloc (lcoefs, res%coefs_multi)
        else
            call move_alloc (lcoefs, res%coefs_multi)
        end if

        if (allocated(res%intercept_multi)) then
            call copy_alloc (lintercept, res%intercept_multi)
        else
            call move_alloc (lintercept, res%intercept_multi)
        end if
    end if

    if (present(ncomp)) ncomp = lncomp

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_pca_1d (conf, X, y, coefs, ncomp, intercept, res, status)
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: coefs
    integer, intent(inout), optional :: ncomp
    real (PREC), intent(out), optional :: intercept
    type (lm_result), intent(inout), optional :: res
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_y
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs
    real (PREC), dimension(:), allocatable :: lintercept

    type (status_t) :: lstatus
    type (lm_config) :: lconf

    lstatus = NF_STATUS_OK

    nullify (ptr_coefs, ptr_y)

    if (present(conf)) lconf = conf
    lconf%trans_y = .false.
    lconf%trans_coefs = .false.

    ptr_y(1:size(y),1:1) => y

    if (present(coefs)) then
        ptr_coefs(1:size(coefs),1:1) => coefs
    end if

    if (present(intercept)) then
        allocate (lintercept(1))
    end if

    call pcr (lconf, X, ptr_y, ptr_coefs, ncomp, lintercept, res, lstatus)

    if (lstatus == NF_STATUS_OK) then
        if (present(res)) then
            if (allocated(res%coefs_multi)) then
                call copy_alloc (res%coefs_multi(:,1), res%coefs)
            end if

            if (allocated(res%intercept_multi)) then
                res%intercept = res%intercept_multi(1)
            end if
        end if

        if (present(intercept)) then
            intercept = lintercept(1)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine



subroutine lm_post_estim (model, X, y, trans_x, rsq, status)
    !*  LM_POST_ESTIM computes common post-estimation statistics such as
    !   R^2 for a given linear model. The model must have been estimated
    !   using before invoking this routine.
    type (lm_result), intent(in) :: model
        !*  Estimated model
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  Array of explanatory (RHS) variables used to estimate model.
        !   Array is assumed to have shape [NOBS, NRHS], unless
        !   TRANS_RHS is present and has value .TRUE.
    real (PREC), intent(in), dimension(:), contiguous :: y
        !*  Vector of observations of the dependent variable used to estimate
        !   model.
    logical, intent(in), optional :: trans_x
    real (PREC), intent(out), optional :: rsq
        !*  If present, contains the model's R^2 on exit.
    type (status_t), intent(out), optional :: status
        !*  Exit code

    logical :: ltrans_x
    real (PREC), dimension(:), allocatable :: resid, y_pred
    character (*), parameter :: NAME = 'POST_ESTIM'

    type (status_t) :: lstatus
    real (PREC) :: var_u, var_y, norm_u
    integer :: Nobs

    lstatus = NF_STATUS_OK

    call set_optional_arg (trans_x, model%conf%trans_x, ltrans_x)

    ! Check whether model has been estimated
    if (.not. allocated(model%coefs)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    call get_regressor_dims (X, trans=ltrans_x, nobs=Nobs)

    call check_cond (size(y) == Nobs, NAME, 'Y: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(rsq)) then

        call compute_resid (lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! Compute variance of residuals
        norm_u = BLAS_NRM2 (model%nobs, resid, 1)
        var_u = norm_u**2.0_PREC / model%nobs

        ! Compute variable of LHS variable Y
        call std (y, s=var_y, dof=0)
        var_y = var_y ** 2.0_PREC

        rsq = 0.0_PREC
        if (var_y > 0.0_PREC) then
            rsq = 1.0_PREC - var_u / var_y
        end if
    end if


100 continue

    if (allocated(y_pred)) deallocate (y_pred)
    if (allocated(resid)) deallocate (resid)

    if (present(status)) status  = lstatus

contains

    subroutine compute_predicted (status)
        type (status_t), intent(out) :: status
        if (.not. allocated (y_pred)) then
            allocate (y_pred(Nobs))
            call predict (X, model%coefs, y_pred, model%intercept, &
                trans_x=ltrans_x, trans_y=.false., trans_coefs=.false., &
                status=status)
        end if
    end subroutine

    subroutine compute_resid (status)
        type (status_t), intent(out) :: status

        call compute_predicted (status)
        if (status /= NF_STATUS_OK) return

        if (.not. allocated (resid)) then
            allocate (resid, source=y)
            resid(:) = resid - y_pred
        end if
    end subroutine

end subroutine



subroutine lm_predict_2d (model, X, Y_pred, trans_x, trans_y, status)
    type (lm_result), intent(in), target :: model
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(out), dimension(:,:), contiguous :: Y_pred
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
    type (status_t), intent(out), optional :: status

    logical :: ltrans_x, ltrans_y
    integer :: Nlhs
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'PREDICT'

    lstatus = NF_STATUS_OK

    call set_optional_arg (trans_x, model%conf%trans_x, ltrans_x)
    call set_optional_arg (trans_y, model%conf%trans_y, ltrans_y)

    call check_cond (allocated(model%coefs_multi), NAME, &
        'MODEL: COEFS_MULTI not allocated', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (allocated(model%intercept_multi), NAME, &
        'MODEL: INTERCEPT_MULTI not allocated', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call get_dims (X, Y_pred, model%coefs_multi, trans_x=ltrans_x, &
        trans_y=ltrans_y, nlhs=Nlhs, status=lstatus)

    call check_cond (model%nlhs == Nlhs, NAME, 'Y_PRED: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call predict (X, model%coefs_multi, Y_pred, intercept=model%intercept_multi, &
        trans_x=ltrans_x, trans_y=ltrans_y, trans_coefs=.false., status=lstatus)

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine lm_predict_1d (model, X, y_pred, trans_x, trans_y, status)
    type (lm_result), intent(in), target :: model
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(out), dimension(:), contiguous :: y_pred
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
        !*  Ignored. Present only for API-compatibility with 2d-variant
        !   of this routine.
    type (status_t), intent(out), optional :: status

    logical :: ltrans_x
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'PREDICT'
    real (PREC), dimension(:), pointer, contiguous :: ptr_coefs

    lstatus = NF_STATUS_OK

    call set_optional_arg (trans_x, model%conf%trans_x, ltrans_x)

    call check_cond (model%nlhs == 1, NAME, 'Y_PRED: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (allocated(model%coefs)) then
        ptr_coefs => model%coefs
    else if (allocated(model%coefs_multi)) then
        ptr_coefs => model%coefs_multi(:,1)
    else
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Prevent intercept being added twice, so pass INTERCEPT argument
    ! only if intercept is not included in COEFS already.
    call predict (X, model%coefs, y_pred, intercept=model%intercept, &
        trans_x=ltrans_x, status=lstatus)

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine predict_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
        trans_y, trans_coefs, status)
    !*  PREDICT computes predicted values for a linear model, ie. it evaluates
    !   some variant of
    !       Y_pred = op(X) * op(B)
    !   or
    !       Y_pred^T = op(B) * op(X)
    !   where all three arrays can be optionally provided in transposed form.
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: coefs
    !*  Coefficient matrix shaped [NCOEFS,NLHS]
    real (PREC), intent(out), dimension(:,:), contiguous :: Y_pred
    real (PREC), intent(in), dimension(:), contiguous, optional :: intercept
    integer, intent(in), dimension(:), contiguous, optional :: irhs
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
    logical, intent(in), optional :: trans_coefs
    type (status_t), intent(out), optional :: status

    logical :: ltrans_y

    call set_optional_arg (trans_y, .false., ltrans_y)

    if (ltrans_y) then
        call predict_trans_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
            trans_coefs, status)
    else
        call predict_notrans_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
            trans_coefs, status)
    end if

end subroutine



subroutine predict_notrans_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
        trans_coefs, status)
    !*  PREDICT_NOTRANS computes predicted values for linear models when
    !   the outcome variables are not transposed, ie. it returns
    !       Y = op(X) * op(B)
    !   where Y has shape [Nobs,Nlhs].
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: coefs
        !*  Coefficient matrix shaped [NCOEFS,NLHS]
    real (PREC), intent(out), dimension(:,:), contiguous :: Y_pred
    real (PREC), intent(in), dimension(:), contiguous, optional :: intercept
    integer, intent(in), dimension(:), contiguous, optional :: irhs
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_coefs
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    logical :: ltrans_x, ltrans_coefs
    integer :: Nrhs, Nlhs, Nobs, Ncoefs, Nconst_coefs, Nrhs_in
    integer :: i, j, jrhs
    character (1) :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real (PREC) :: beta, b0, b_ij
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_b
    real (PREC), dimension(:,:), allocatable :: X_in
    logical :: deallocate_coefs

    character (*), parameter :: NAME = 'PREDICT'

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (trans_x, .false., ltrans_x)
    call set_optional_arg (trans_coefs, .false., ltrans_coefs)

    call get_dims (X, Y_pred, coefs, trans_x=ltrans_x, &
        trans_y=.false., trans_coefs=ltrans_coefs, nrhs=Nrhs, nlhs=Nlhs, &
        nobs=Nobs, ncoefs=Ncoefs, nconst_coefs=Nconst_coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(intercept)) then
        call check_cond (size(intercept) == Nlhs, NAME, &
            'INTERCEPT: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Add intercept ---

    if (present(intercept)) then
        do j = 1, Nlhs
            Y_pred(:,j) = intercept(j)
        end do
        ! Add X'B without intercept in GEMM, GEMV call
        beta = 1.0
    else if (Nconst_coefs == 1) then
        ! Initialize Y with intercepts from first row in COEFS
        do j = 1, Nlhs
            if (ltrans_coefs) then
                b0 = coefs(j,1)
            else
                b0 = coefs(1,j)
            end if
            Y_pred(:,j) = b0
        end do

        ! Add X'B without intercept in GEMM calls
        beta = 1.0
    else
        beta = 0.0
    end if

    ! --- Add XB in some form ---

    if (present(irhs)) then
        Nrhs_in = size(irhs)
    else
        Nrhs_in = Nrhs
    end if

    ! Exit immediately for degenerate matrix X
    if (Nrhs_in == 0) goto 100

    ! We need to account for the following variants:
    !   1.  X and COEFS not transposed: Y += X*B(i:,:)
    !   2.  X transposed, COEFS not transposed: Y += X^T * B(i:,:)
    !   3.  X not transposed, COEFS transposed: Y += X * B(:,i:)^T
    !   4.  X transposed, COEFS transposed: Y += X^T * B(:,i:)^T
    !
    !   where the initial row i depends on whether COEFS includes
    !   an intercept that is not in X.

    if (present(irhs)) then
        ! Include only selected RHS indices.

        ! If there is no intercept, Y_PREC needs to be initialized to zero.
        if (.not. (Nconst_coefs == 1 .or. present(intercept))) then
            Y_pred(:,:) = 0.0
        end if

        if (.not. ltrans_x) then
            ! X not tranposed: simple case with mostly memory-contiguous access,
            ! so that we can just loop through columns of Y and X.
            do i = 1, Nlhs
                do j = 1, Nrhs_in
                    jrhs = irhs(j)
                    if (ltrans_coefs) then
                        b_ij = coefs(i,Nconst_coefs+jrhs)
                    else
                        b_ij = coefs(Nconst_coefs+jrhs,i)
                    end if

                    if (b_ij == 0.0_PREC) cycle

                    call BLAS_AXPY (Nobs, b_ij, X(:,jrhs), 1, Y_pred(:,i), 1)
                end do
            end do
        else
            ! We need to compute
            !   Y = X(irhs,:)^T * op(B)(Nconst+irhs,:)
            ! which cannot be accomplished in a contiguous fashion in any
            ! reasonable way. Instead create a contiguous versions of X.
            ! and compute
            !   Y = X_in * op(B)(Nconst+irhs,:)

            ! Allocate packed X and COEFS
            allocate (X_in(Nobs, Nrhs_in))

            call replay_transform (X, X_in, trans=.true., ikeep=irhs)

            do i = 1, Nlhs
                do j = 1, Nrhs_in
                    jrhs = irhs(j)
                    if (ltrans_coefs) then
                        b_ij = coefs(i,Nconst_coefs+jrhs)
                    else
                        b_ij = coefs(Nconst_coefs+jrhs,i)
                    end if

                    if (b_ij == 0.0_PREC) cycle

                    call BLAS_AXPY (Nobs, b_ij, X_in(:,j), 1, Y_pred(:,i), 1)
                end do
            end do
        end if
    else
        ! Include all RHS variables.


        ! Y is not transposed, we compute some variant of
        ! Y = op(X) * op(B)(1+Nconst:,:)
        ldc = Nobs
        m = Nobs
        n = Nlhs
        k = Nrhs

        if (ltrans_x) then
            ! Case 2 + 4
            ! X provided in transposed form, needs to be transposed back to
            ! [Nobs,Nrhs]
            transa = 'T'
            lda = Nrhs
        else
            ! X provided in non-transposed form, remains as is
            ! Case 1 + 3
            transa = 'N'
            lda = Nobs
        end if

        ! Shift initial relevant row or column depending on whether intecept is
        ! (additionally) stored in COEFS.
        deallocate_coefs = .false.
        if (ltrans_coefs) then
            ! Case  3 + 4
            ! COEFS transposed by caller, needs to be transposed back
            ! to [NCOEFS,NLHS]
            transb = 'T'
            ldb = Nlhs
            ptr_b => coefs(:,1+Nconst_coefs:Ncoefs)
        else
            ! Case 1 + 2
            ! COEFS provided in non-transposed form, remains as is
            transb = 'N'
            ldb = Nrhs
            if (Nconst_coefs == 1) then
                deallocate_coefs = .true.
                allocate (ptr_b(Nrhs,Nlhs))
                do i = 1, Nlhs
                    call BLAS_COPY (Nrhs, coefs(2:Ncoefs,i), 1, ptr_b(:,i), 1)
                end do
            else
                ptr_b => coefs
            end if
        end if

        call BLAS_GEMM (transa, transb, m, n, k, 1.0_PREC, X, lda, &
            ptr_b, ldb, beta, Y_pred, ldc)

        if (deallocate_coefs) then
            deallocate (ptr_b)
        end if
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine predict_trans_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
        trans_coefs, status)
    !*  PREDICT_TRANS implements prediction for linear models when
    !   the outcome array is transposed, ie.
    !       Y^T = op(B) * op(X)
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: coefs
    !*  Coefficient matrix shaped [NCOEFS,NLHS]
    real (PREC), intent(out), dimension(:,:), contiguous :: Y_pred
    real (PREC), intent(in), dimension(:), contiguous, optional :: intercept
    integer, intent(in), dimension(:), contiguous, optional :: irhs
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_coefs
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    logical :: ltrans_x, ltrans_coefs
    integer :: Nrhs, Nlhs, Nobs, Ncoefs, Nconst_coefs, Nrhs_in
    integer :: i, j, jrhs
    character (1) :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real (PREC) :: beta
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_b
    real (PREC), dimension(:,:), allocatable :: X_in_T, coefs_in_T
    character (*), parameter :: NAME = 'PREDICT'
    logical :: deallocate_coefs

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (trans_x, .false., ltrans_x)
    call set_optional_arg (trans_coefs, .false., ltrans_coefs)

    call get_dims (X, Y_pred, coefs, trans_x=ltrans_x, trans_y=.true., &
        trans_coefs=ltrans_coefs, nrhs=Nrhs, nlhs=Nlhs, nobs=Nobs, &
        ncoefs=Ncoefs, nconst_coefs=Nconst_coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(intercept)) then
        call check_cond (size(intercept) == Nlhs, NAME, &
            'INTERCEPT: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Add intercept ---

    if (present(intercept)) then
        Y_pred(:,:) = spread(intercept, dim=2, ncopies=Nobs)
        ! Add X'B without intercept in GEMM, GEMV call
        beta = 1.0
    else if (Nconst_coefs == 1) then
        ! Initialize Y with intercepts from first row in COEFS
        if (ltrans_coefs) then
            Y_pred(:,:) = spread(coefs(:,1), dim=2, ncopies=Nobs)
        else
            Y_pred(:,:) = spread(coefs(1,:), dim=2, ncopies=Nobs)
        end if

        ! Add X'B without intercept in GEMM calls
        beta = 1.0
    else
        beta = 0.0
    end if

    ! --- Add XB in some form ---

    if (present(irhs)) then
        Nrhs_in = size(irhs)
    else
        Nrhs_in = Nrhs
    end if

    ! Exit immediately for degenerate matrix X
    if (Nrhs_in == 0) goto 100

    ! We need to account for the following variants:
    !   1.  X and COEFS not transposed: Y^T += B(i:,:)^T * X^T
    !   2.  X transposed, COEFS not transposed: Y^T += B(i:,:)^T * X
    !   3.  X not transposed, COEFS transposed: Y^T += B(:,i:) * X^T
    !   4.  X transposed, COEFS transposed: Y^T += B(:,i:) * X
    !
    !   where the initial row i depends on whether COEFS includes
    !   an intercept that is not in X.

    if (present(irhs)) then
        ! Include only selected RHS indices.
        ! This is impossible to solve in a way that uses contiguous
        ! memory access, so we first create contiguous packed versions
        ! of X and COEFS.
        !
        ! We need to compute Y^T += B^T * X^T
        ! so we need to create transposed arrays of B and X, if these
        ! are not already on entry.

        ! Allocate packed X and COEFS
        allocate (X_in_T(Nrhs_in, Nobs))
        allocate (coefs_in_T(Nlhs,Nrhs_in))

        if (ltrans_x) then
            do i = 1, Nobs
                do j = 1, Nrhs_in
                    jrhs = irhs(j)
                    X_in_T(j,i) = X(jrhs,i)
                end do
            end do
        else
            do i = 1, Nobs
                do j = 1, Nrhs_in
                    jrhs = irhs(j)
                    X_in_T(j,i) = X(i,jrhs)
                end do
            end do
        end if

        if (ltrans_coefs) then
            ! Transposed COEFS with shape [NLHS,NCOEFS]
            do j = 1, Nrhs_in
                jrhs = irhs(j)
                coefs_in_T(:,j) = coefs(1:Nlhs,Nconst_coefs+jrhs)
            end do
        else
            ! Non-transposed COEFS with shape [NCOEFS,NLHS]
            do j = 1, Nrhs_in
                jrhs = irhs(j)
                coefs_in_T(:,j) = coefs(Nconst_coefs+jrhs,1:Nlhs)
            end do
        end if

        call BLAS_GEMM ('N', 'N', Nlhs, Nobs, Nrhs_in, 1.0_PREC, coefs_in_T, &
            Nlhs, X_in_T, Nrhs_in, beta, Y_pred, Nlhs)
    else
        ! Include all RHS variables.

        ldc = Nlhs
        m = Nlhs
        n = Nobs
        k = Nrhs

        if (ltrans_x) then
            ! Case 2 + 4
            ! X provided in transposed form, remains as is
            transb = 'N'
            ldb = Nrhs
        else
            ! Case 1 + 3
            ! X provided in non-transposed form, needs to be transposed to
            ! shape [Nrhs,Nobs]
            transb = 'T'
            ldb = Nobs
        end if

        deallocate_coefs = .false.
        if (ltrans_coefs) then
            ! Case 3 + 4
            ! Already transposed by caller, remains as is
            transa = 'N'
            lda = Nlhs
            ptr_b => coefs(:,1+Nconst_coefs:Ncoefs)
        else
            ! Case 1 + 2
            ! COEFS provided in non-transposed form, needs to be transposed
            ! to shape [Nlhs,Ncoefs]
            transa = 'T'
            lda = Nrhs
            if (Nconst_coefs == 1) then
                deallocate_coefs = .true.
                allocate (ptr_b(Nrhs,Nlhs))
                do i = 1, Nlhs
                    call BLAS_COPY (Nrhs, coefs(2:Ncoefs,i), 1, ptr_b(:,i), 1)
                end do
            else
                ptr_b => coefs
            end if
        end if

        call BLAS_GEMM (transa, transb, m, n, k, 1.0_PREC, ptr_b, lda, &
            X, ldb, beta, Y_pred, ldc)

        if (deallocate_coefs) then
            deallocate (ptr_b)
        end if
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine predict_1d (X, coefs, y_pred, intercept, irhs, trans_x, trans_y, &
        trans_coefs, status)
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: coefs
        !*  Coefficient vector of length NCOEFS, where NCOEFS = NRHS + 1
        !   if ADD_INTERCEPT=.TRUE., and NCOEFS = NRHS otherwise.
    real (PREC), intent(out), dimension(:), contiguous, target :: y_pred
    real (PREC), intent(in), optional :: intercept
    integer, intent(in), dimension(:), contiguous, optional :: irhs
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
        !*  Ignored. Present only for API compatibility with 2d variant of this
        !   routine.
    logical, intent(in), optional :: trans_coefs
        !*  Ignored. Present only for API compatibility with 2d variant of this
        !   routine.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    logical :: ltrans_x, ltrans_y, ltrans_coefs
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_y
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs
    real (PREC), dimension(:), allocatable :: lintercept

    call set_optional_arg (trans_x, .false., ltrans_x)

    if (ltrans_x) then
        ! The optimal path to take is to compute
        !   Y^T = B^T * X^T
        ! at least if IRHS is not present, since that is only one GEMM call
        ! and no additional memory allocations.
        ptr_coefs(1:1,1:size(coefs)) => coefs
        ptr_y(1:1,1:size(y_pred)) => y_pred
        ltrans_y = .true.
        ltrans_coefs = .true.
    else
        ! The optimal path to take is
        !   Y = X * B
        ! since that causes no additional allocations and only uses contiguous
        ! memory access.
        ptr_coefs(1:size(coefs),1:1) => coefs
        ptr_y(1:size(y_pred),1:1) => y_pred
        ltrans_y = .false.
        ltrans_coefs = .false.
    end if

    if (present(intercept)) then
        allocate (lintercept(1), source=intercept)
    end if

    call predict (X, ptr_coefs, ptr_y, intercept=lintercept, irhs=irhs, &
        trans_x=trans_x, trans_y=ltrans_y, trans_coefs=ltrans_coefs, status=status)

end subroutine
