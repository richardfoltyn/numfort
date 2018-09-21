


subroutine __APPEND(root_broyden_check_input,__PREC) (maxiter, maxfev, &
        tol, xtol, rstep, xstep, status, msg)
    integer, parameter :: PREC = __PREC
    integer, intent(in):: maxiter, maxfev
    real (PREC), intent(in) :: tol, xtol
    real (PREC), intent(in), optional :: rstep, xstep
    type (status_t), intent(inout) :: status
    character (*), intent(inout) :: msg

    status = NF_STATUS_OK

    call check_nonneg (0, maxiter, "maxiter", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_nonneg (0, maxfev, "maxfev", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, tol, "tol", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, xtol, "xtol", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, rstep, "rstep", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, xstep, "xstep", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG
end subroutine


recursive subroutine __APPEND(root_broyden,__PREC) (fcn, x, ndiff, tol, xtol, &
        maxiter, maxfev, rstep, xstep, dstep, iprint, work, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fvv_fcn,__PREC)) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    real (PREC), intent(in), optional :: dstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
        !*  If present, debug info will be printed according to the value
        !   of IPRINT.
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_vv,__PREC)) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn=fcn, eps=dstep)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, res)

end subroutine


recursive subroutine __APPEND(root_broyden_jac,__PREC) (fcn, fjac, x, tol, &
        xtol, maxiter, maxfev, rstep, xstep, iprint, work, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fvv_fcn,__PREC)) :: fcn
    procedure (__APPEND(fvv_jac,__PREC)) :: fjac
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_vv,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn, jac=fjac)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, res)

end subroutine

recursive subroutine __APPEND(root_broyden_fcn_jac_opt,__PREC) (fcn, x, tol, &
        xtol, maxiter, maxfev, rstep, xstep, iprint, work, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fvv_fcn_jac_opt,__PREC)) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_vv,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac_opt=fcn)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, res)

end subroutine


recursive subroutine __APPEND(root_broyden_args,__PREC) (fcn, x, args, ndiff, &
        tol, xtol, maxiter, maxfev, rstep, xstep, dstep, iprint, work, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fvv_fcn_args,__PREC)) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    real (PREC), intent(in), optional :: dstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_vv,__PREC)) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn_args=fcn, args=args, eps=dstep)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, res)

end subroutine


recursive subroutine __APPEND(root_broyden_jac_args,__PREC) (fcn, fjac, x, &
        args, tol, xtol, maxiter, maxfev, rstep, xstep, iprint, work, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fvv_fcn_args,__PREC)) :: fcn
    procedure (__APPEND(fvv_jac_args,__PREC)) :: fjac
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_vv,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_args=fcn, jac_args=fjac, args=args)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, res)

end subroutine



recursive subroutine __APPEND(root_broyden_fcn_jac_opt_args,__PREC) (fcn, x, &
        args, tol, xtol, maxiter, maxfev, rstep, xstep, iprint, work, res)
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fvv_fcn_jac_opt_args,__PREC)) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_vv,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac_opt_args=fcn, args=args)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, res)

end subroutine


recursive subroutine __APPEND(root_broyden_impl,__PREC) (fcn, x, tol, xtol, &
        maxiter, maxfev, rstep, xstep, iprint, work, res)

    integer, parameter :: PREC = __PREC

    type (__APPEND(fwrapper_vv,__PREC)), intent(inout) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
        !*  Max. number of function evaluations (includes evaluations of
        !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
        !*  Max. step size in search direction, relative to the current
        !   point (default: unbounded)
    real (PREC), intent(in), optional :: xstep
        !*  Max. absolute step size in search direction (default: unbounded)
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    type (__APPEND(workspace,__PREC)), intent(inout), optional, target :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), optional, target :: res

    real (PREC) :: ltol, lxtol
    real (PREC) :: dx_scale, denom, nrm, nrmp1, nrm_last, nrm_upd
    integer :: lmaxiter, lmaxfev, k, n, i, nrwrk, niwrk, liprint
    integer :: lwork_inv, liwork_inv
        !   Workspace array sizes use for INV routine
    integer, dimension(2) :: shp2d
    real (PREC), dimension(:), pointer, contiguous :: fx, fxlast, dx, dfx
    real (PREC), dimension(:), pointer, contiguous :: vec1, vec2
        !   Temporary arrays used for various BLAS routines
    real (PREC), dimension(:), pointer, contiguous :: x_ls, fx_ls
        !   Working arrays used for line search
    real (PREC), dimension(:,:), pointer, contiguous :: jac, jac_inv
    real (PREC), dimension(:), pointer, contiguous :: rwork_inv
    integer, dimension(:), pointer, contiguous :: iwork_inv

    type (__APPEND(workspace,__PREC)), pointer :: ptr_work
    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res
    type (status_t) :: status

    ! Arguments for GEMV
    character (1) :: trans
    real (PREC) :: alpha, beta
    integer, parameter :: incx = 1, incy = 1
    integer :: lda, m
    
    status = NF_STATUS_OK
    nullify (ptr_work, ptr_res)

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    ! Default arguments
    lmaxiter = 100
    ltol = sqrt(epsilon(0.0_PREC))
    lxtol = sqrt(epsilon(0.0_PREC))
    liprint = PRINT_NONE

    ! Overwrite defaults with optionally provided user arguments
    ! Note: RSTEP and XSTEP have no default values, we skip limiting the
    ! step size completly whenever they are not present.
    if (present(maxiter)) lmaxiter = maxiter
    if (present(tol)) ltol = tol
    if (present(xtol)) lxtol = xtol
    if (present(iprint)) liprint = iprint

    ! Max. number of function evaluations given by max. iteration count and
    ! max. number of backtracking steps performed during linesearch.
    ! Add 1 for initial function evaluation, 1 to ensure that root finder
    ! will not exit due to NFEV exceeding MAXFUN when it in fact should
    ! terminate due to exceeding MAXITER.
    lmaxfev = lmaxiter * LINESEARCH_MAX_STEPS + 2
    if (fcn%num_diff) then
        ! If numerical differentiation is performed
        lmaxfev = lmaxfev + size(x)
    end if
    if (present(maxfev)) lmaxfev = maxfev

    ! Validate inputs
    call root_broyden_check_input (lmaxiter, lmaxfev, ltol, lxtol, &
        rstep, xstep, status, ptr_res%msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    n = size(x)
    k = 0

    ! Arguments to GEMV that never change
    m = n
    lda = n

    ! Workspace array size
    nrwrk = 2*n*n + 8*n
    ! Working arrays needed for INV
    call inv_work_query (n, lwork_inv, liwork_inv)
    nrwrk = nrwrk + lwork_inv
    niwrk = liwork_inv

    ! Allocate workspace arrays
    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal state in workspace object, in particular index offsets
    ! (this does not deallocate working arrays)
    call workspace_reset (ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk, niwrk=niwrk)

    call workspace_get_ptr (ptr_work, n, vec1)
    call workspace_get_ptr (ptr_work, n, vec2)
    call workspace_get_ptr (ptr_work, n, fx)
    call workspace_get_ptr (ptr_work, n, fxlast)
    call workspace_get_ptr (ptr_work, n, dx)
    call workspace_get_ptr (ptr_work, n, dfx)
    call workspace_get_ptr (ptr_work, n, fx_ls)
    call workspace_get_ptr (ptr_work, n, x_ls)

    shp2d = n
    call workspace_get_ptr (ptr_work, shp2d, jac)

    ! Working arrays for INV
    call workspace_get_ptr (ptr_work, shp2d, jac_inv)
    call workspace_get_ptr (ptr_work, n, rwork_inv)
    call workspace_get_ptr (ptr_work, n, iwork_inv)

    if (liprint /= PRINT_NONE) then
        print '(tr1, a)', "ROOT_BROYDEN: routine start"
    end if

    call dispatch (fcn, x, fxlast)

    if (.not. all(ieee_is_finite(fxlast))) then
        ptr_res%msg = 'Invalid function value encountered'
        status = NF_STATUS_INVALID_STATE
        goto 100
    end if

    ! Check whether initial point satisfied convergence criterion
    nrm_last = norm2(fxlast)

    if (liprint /= PRINT_NONE) then
        print '(tr1, a)', "ROOT_BROYDEN: initial values"
        print 200, "x: ", x
        print 200, "f(x): ", fxlast
        print 201, "2-norm: ", nrm_last
    end if

    if (nrm_last < ltol) then
        ptr_res%msg = 'Convergence achieved, func. value smaller than tol'
        status = NF_STATUS_OK
        goto 100
    end if

    ! Compute Jacobian at initial point. Reuse f(X) evaluated above in
    ! case of numerical differentiation.
    call dispatch_jac (fcn, x, jac, fxlast)

    if (.not. all(ieee_is_finite(jac))) then
        ptr_res%msg = 'Invalid value in initial Jacobian encountered'
        status = NF_STATUS_INVALID_STATE
        goto 100
    end if

    if (iand(liprint, PRINT_JAC) == PRINT_JAC) then
        print '(tr1, a)', "ROOT_BROYDEN: initial Jacobian"
        do i = 1, n
            print 200, "", jac(i,:)
        end do
    end if

    ! Initial inverted Jacobian
    call inv (jac, jac_inv, rwork=rwork_inv, iwork=iwork_inv, status=status)
    if (status /= NF_STATUS_OK) then
        ptr_res%msg = "Could not invert initial Jacobian"
        goto 100
    end if

   do k = 1, lmaxiter

        ! 1. Find next candidate point
        alpha = -1.0_PREC
        beta = 0.0_PREC
        trans = 'N'
        ! Compute dx_k = - J^{-1}_{k-1} * f(x_{k-1})
        call GEMV (trans, m, n, alpha, jac_inv, lda, fxlast, incx, beta, dx, incy)

        dx_scale = 1.0_PREC
        if (present(rstep)) then
            ! Perform "dampened" update step, allowing point-wise steps of
            ! at most rstep * abs(x_k).
            ! Standard Broyden would set x_k = x_{k-1} + dx_k
            do i = 1, n
                dx_scale = min(dx_scale, rstep * abs(x(i)) / abs(dx(i)))
            end do
        end if

        if (present(xstep)) then
            ! Check that step size does not violate max. absolute step
            ! size for any dimension.
            do i = 1, n
                dx_scale = min(dx_scale, xstep / abs(dx(i)))
            end do
        end if

        if (iand(liprint, PRINT_STEP) == PRINT_STEP) then
            print '(tr1,a)', "ROOT_BROYDEN: STEP direction"
            print 200, "Unrestricted step direction: ", dx
            print 201, "Rescaling factor: ", dx_scale
        end if

        ! Rescale step size such that change in x is at most EPSX percent
        ! for each element.
        dx(:) = dx_scale * dx

        call dumb_line_search (fcn, x, nrm_last, liprint, x_ls, fx_ls, dx, fx)

        if (.not. all(ieee_is_finite(fx))) then
            ptr_res%msg = 'Invalid function value encountered during line search'
            status = NF_STATUS_INVALID_STATE
            goto 100
        end if

         ! 2. Check convergence in terms of function values
        nrm_upd = norm2(fx)
        if (nrm_upd < ltol) then
            ptr_res%msg = 'Convergence achieved, func. value smaller than tol'
            status = NF_STATUS_OK
            ! Code below expects the final function value to be stored in FXLAST
            call DCOPY (n, fx, 1, fxlast, 1)
            x(:) = x + dx

            goto 100
        end if

        ! 3. Check if x_k and x_{k-1} are sufficiently close to each other to
        ! terminate algorithm
        nrm = norm2(dx)
        nrmp1 = norm2(x) + 1.0_PREC
        if (nrm < lxtol * nrmp1) then
            ptr_res%msg = 'Not making sufficient progress, change in x smaller than xtol'
            status = NF_STATUS_NOT_CONVERGED

            ! Determine whether to return (X,FXLAST) or (X+DX,FX) is a
            ! better approximation
            if (nrm_last > nrm_upd) then
                x(:) = x + dx
                call DCOPY (n, fx, 1, fxlast, 1)
            end if

            goto 100
        end if

        ! 4. Check whether max. number of func. evaluations was exceeded.
        if (fcn%nfev >= lmaxfev) then
            ! Set corresponding exit status.
            ! At this point the last "best" guess for the root is stored in X,
            ! the corresponding function value in FXLAST.
            ptr_res%msg = 'Max. number of function evaluations exceeded'
            status = NF_STATUS_MAX_EVAL

            ! Determine whether to return (X,FXLAST) or (X+DX,FX) is a
            ! better approximation
            if (nrm_last > nrm_upd) then
                x(:) = x + dx
                call DCOPY (n, fx, 1, fxlast, 1)
            end if

            goto 100
        end if

        ! Update X with the "best" DX identified by line search
        x(:) = x + dx

        ! Compute difference to last objective (fx contains the best linesearch
        ! result)
        dfx(:) = fx - fxlast

        ! 4. Update inverse of J_k
        ! Compute J^{-1}_{k-1} dfx_k
        alpha = 1.0_PREC
        beta = 0.0_PREC
        trans = 'N'
        call GEMV (trans, m, n, alpha, jac_inv, lda, dfx, incx, beta, vec1, incy)
        ! Denominator in expression for J^{-1}_k, dx_k^T J^{-1}_{k-1} dfx_k
        denom = DOT (n, vec1, incx, dx, incy)

        ! Compute transpose of dx_k^T * J^{-1}_{k-1}
        trans = 'T'
        alpha = 1.0_PREC
        beta = 0.0_PREC
        call GEMV (trans, m, n, alpha, jac_inv, lda, dx, incx, beta, vec2, incy)

        ! Compute X argument in call to DGER, ie. dx_k - J^{-1}_{k-1} dfx_k
        ! Set VEC1 = - J^{-1}_{k-1} dfx_k
        alpha = -1.0_PREC
        call SCAL (n, alpha, vec1, incx)
        ! Add dx_k to VEC1
        alpha = 1.0_PREC
        call AXPY (n, alpha, dx, incx, vec1, incy)
        ! Compute update to inverse Jacobian using BLAS rank-1 update routine
        ! which adds alpha * VEC1 * VEC2^T to J^{-1}_{k-1}
        alpha = 1.0_PREC / denom
        call GER (m, n, alpha, vec1, incx, vec2, incy, jac_inv, lda)

        ! Update for next iteration
        call DCOPY (n, fx, 1, fxlast, 1)
        nrm_last = nrm_upd
    end do

    ptr_res%msg = 'Max. number of iterations exceeded'
    status = NF_STATUS_MAX_ITER
    ! Correct number of iterations
    k = lmaxiter

100 continue

    if (present(res)) then
        if (associated(fxlast)) then
            call result_update (ptr_res, x, fxlast, status, nit=k, nfev=fcn%nfev)
        else
            ! FXLAST may not be associated if routine exists prematurely
            ! due to errors
            call result_update (ptr_res, x, status=status, nit=k, nfev=fcn%nfev)
        end if
    end if

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

    if (liprint /= PRINT_NONE) then
        print '(tr1, a)', "ROOT_BROYDEN: routine complete"
    end if

200 format (tr3, a, *(es12.5e2,:,", "))
201 format (tr3, a, es12.5e2)

end subroutine



recursive subroutine __APPEND(dumb_line_search,__PREC) (fcn, x, nrm, iprint, &
        x_ls, fx_ls, dx, fx)
    !*  DUMB_LINE_SEARCH uses backtracking to identify the step size in a
    !   given direction that yields a lower Euclidean norm of the objective
    !   function than at a given point X.
    !
    !   If the Euclidean norm is larger at all candidate points in direction
    !   DX, the least bad step size is returned.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vv,__PREC)), intent(inout) :: fcn
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in) :: nrm
        !*  Euclidean norm of the last function value, |FX|_2
    integer (NF_ENUM_KIND), intent(in) :: iprint
    real (PREC), intent(inout), dimension(:), contiguous :: x_ls, fx_ls
        !*  Working arrays to store X, FX evaluated during line search
    real (PREC), intent(inout), dimension(:) :: dx
        !*  On entry, contains the max. step size in a desired direction.
        !   On exit, the value is updated to contain the "optimal" step size.
    real (PREC), intent(inout), dimension(:), contiguous :: fx
        !*  On entry, contains the function value evaluated at X. On exit,
        !   contains the function value at X + DX where DX is the updated
        !   value computed by the routine.

    integer :: i, n
    real (PREC) :: step_size, nrm_ls, nrm_ls_best, step_best

    n = size(x)

    ! Line search
    nrm_ls_best = huge(1.0_PREC)
    step_best = 0.0_PREC

    if (iand(iprint, PRINT_LSEARCH) == PRINT_LSEARCH) then
        print '(tr1, a)', "ROOT_BROYDEN: LINE SEARCH routine"
    end if

    do i = 1, LINESEARCH_MAX_STEPS
        step_size = (LINESEARCH_MAX_STEPS-i+1.0_PREC)/LINESEARCH_MAX_STEPS
        x_ls(:) = x + step_size * dx

        call dispatch (fcn, x_ls, fx_ls)
        nrm_ls = norm2(fx_ls)

        if (iand(iprint, PRINT_LSEARCH) == PRINT_LSEARCH) then
            if (i > 1) then
                print '(tr3, a, i0)', "ROOT_BROYDEN: BACKTRACKING step #", i-1
            end if
            print 200, "x: ", x_ls
            print 200, "f(x): ", fx_ls
            print 201, "2-norm: ", nrm_ls
        end if

        if (nrm_ls <= nrm) then
            fx(:) = fx_ls
            dx = dx * step_size
            goto 100
        else if (nrm_ls < nrm_ls_best) then
            fx(:) = fx_ls
            step_best = step_size
            nrm_ls_best = nrm_ls
        end if
    end do

    ! No step size with decreasing objective found, take the least bad
    dx = step_best * dx

100 continue

    if (iand(iprint, PRINT_LSEARCH) == PRINT_LSEARCH) then
        print '(tr1, a)', "ROOT_BROYDEN: LINE SEARCH complete"
    end if

200 format (tr3, a, *(es12.5e2, :, ", "))
201 format (tr3, a, es12.5e2)

end subroutine

