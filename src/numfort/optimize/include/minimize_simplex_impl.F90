
subroutine __APPEND(check_input,__PREC) (tol, maxfun, status, msg)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: tol
    integer, intent(in) :: maxfun
    type (status_t), intent(inout) :: status
    character (*), intent(inout) :: msg

    status = NF_STATUS_OK

    call check_positive (0, maxfun, "maxfun", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, tol, "tol", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



subroutine __APPEND(minimize_simplex,__PREC) (fcn, x, tol, maxfun, quad, &
        iprint, work, res)

    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fvs_fcn,__PREC)) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxfun
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    logical, intent(in), optional :: quad
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_vs,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn)

    call simplex_impl (fwrapper, x, tol, maxfun, quad, iprint, work, res)
end subroutine


subroutine __APPEND(minimize_simplex_args,__PREC) (fcn, x, args, tol, maxfun, &
        quad, iprint, work, res)

    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fvs_fcn_args,__PREC)) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxfun
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    logical, intent(in), optional :: quad
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), optional :: res

    type (__APPEND(fwrapper_vs,__PREC)) :: fwrapper

    call wrap_procedure (fwrapper, fcn_args=fcn, args=args)

    call simplex_impl (fwrapper, x, tol, maxfun, quad, iprint, work, res)
end subroutine


subroutine __APPEND(simplex_impl,__PREC) (fcn, x, tol, maxfun, quad, &
        iprint, work, res)

    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vs,__PREC)), intent(inout) :: fcn

    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxfun
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    logical, intent(in), optional :: quad
    type (__APPEND(workspace,__PREC)), intent(inout), target, optional :: work
    type (__APPEND(optim_result,__PREC)), intent(inout), target, optional :: res

    integer :: n, lmaxfun, nloop, iquad, liprint, ifault, nrwrk
    logical :: lquad
    real (PREC) :: delta_nonzero, delta_zero, simp, fopt, ltol
    real (PREC), dimension(:), pointer, contiguous :: ptr_step, ptr_var

    type (__APPEND(workspace,__PREC)), pointer :: ptr_work
    type (__APPEND(optim_result,__PREC)), pointer :: ptr_res

    nullify (ptr_work, ptr_res)
    nullify (ptr_step, ptr_var)

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    n = size(x)
    ltol = 1.0e-4_PREC
    lmaxfun = 100 * n
    lquad = .true.
    iquad = 0

    if (present(tol)) ltol = tol
    if (present(maxfun)) lmaxfun = maxfun
    if (present(quad)) lquad = quad
    if (lquad) iquad = 1

    ! disable diagnostic printing
    liprint = map_iprint(NF_PRINT_NONE)
    if (present(iprint)) liprint = map_iprint (iprint)

    call check_input (ltol, lmaxfun, ptr_res%status, ptr_res%msg)
    if (NF_STATUS_INVALID_ARG .in. ptr_res%status) goto 100

    nrwrk = 2 * n

    ! Allocate workspace arrays
    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal state in workspace object, in particular index offsets
    ! (this does not deallocate working arrays)
    call workspace_reset (ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk)

    call workspace_get_ptr (ptr_work, n, ptr_step)
    call workspace_get_ptr (ptr_work, n, ptr_var)

    ! set up initial step size
    ! take some sensible values as used in scipy implementation
    delta_zero = 0.00025_PREC
    delta_nonzero = 0.05_PREC
    where (x == 0.0_PREC)
        ptr_step = delta_zero
    else where
        ptr_step = x * delta_nonzero
    end where

    nloop = 2 * n
    simp = 1.0e-6_PREC

    call minim (x, ptr_step, n, fopt, lmaxfun, liprint, ltol, nloop, iquad, &
        simp, ptr_var, fcn_wrapper, ifault)

    ! map mimum status into numfort_optimize status
    call map_ifault (ifault, ptr_res%status, ptr_res%msg)

100 continue

    if (present(res)) then
        call result_update (ptr_res, x, fopt, nfev=fcn%nfev)
    end if

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

    contains

    subroutine fcn_wrapper (x, fx)
        real (PREC), intent(in), dimension(:), contiguous :: x
        real (PREC), intent(out) :: fx

        call dispatch (fcn, x, fx)
    end subroutine

end subroutine
