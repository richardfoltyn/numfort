

subroutine __APPEND(polyfit_check_input,__PREC) (x, y, deg, coefs, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    integer, intent(in) :: deg
    real (PREC), intent(in), dimension(:,:) :: coefs
    type (status_t), intent(out) :: status

    integer :: n

    n = size(x)
    status = NF_STATUS_INVALID_ARG

    if (n /= size(y,1)) return
    ! Cannot fit data to polynomial if there are not enough data points for
    ! requested polynomial degree.
    if (deg >= n) return
    if (deg < 0) return

    if (size(coefs,1) < (deg + 1)) return
    if (size(coefs,2) < size(y,2)) return

    status = NF_STATUS_OK
end subroutine


subroutine __APPEND(polyfit_1d,__PREC) (x, y, deg, coefs, work, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:), target :: y
    integer, intent(in) :: deg
    real (PREC), intent(out), dimension(:), contiguous, target :: coefs
    type (__APPEND(workspace,__PREC)), optional :: work
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer :: ptr_y, ptr_coefs

    ptr_y(1:size(y),1:1) => y
    ptr_coefs(1:size(coefs),1:1) => coefs

    call polyfit (x, ptr_y, deg, ptr_coefs, work, status)
end subroutine


subroutine __APPEND(polyfit,__PREC) (x, y, deg, coefs, work, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    integer, intent(in) :: deg
    real (PREC), intent(out), dimension(:,:), contiguous :: coefs
    type (__APPEND(workspace,__PREC)), optional :: work
    type (status_t), intent(out), optional :: status

    type (__APPEND(workspace,__PREC)), pointer :: ptr_work
    type (status_t) :: lstatus

    integer :: n

    lstatus = NF_STATUS_OK

    call polyfit_check_input (x, y, deg, coefs, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    call assert_alloc_ptr (work, ptr_work)

    n = size(x)

    if (deg == (n-1)) then
        call polyfit_exact (x, y, deg, work, coefs, lstatus)
    else
        lstatus = NF_STATUS_UNSUPPORTED_OP
        goto 100
    end if

100 continue

    if (present(status)) status = lstatus

    call assert_dealloc_ptr (work, ptr_work)

end subroutine



subroutine __APPEND(polyfit_exact,__PREC) (x, y, deg, work, coefs, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    integer, intent(in) :: deg
    type (__APPEND(workspace,__PREC)) :: work
    real (PREC), intent(out), dimension(:,:), contiguous :: coefs
    type (status_t), intent(out) :: status

    ! Arguments for calling GESV
    integer :: n, info, nrhs, lda, ldb
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_a
    integer, dimension(:), pointer :: ptr_ipiv

    integer, dimension(2) :: shp

    n = size(x)
    nrhs = size(y, 2)
    lda = n
    ldb = deg + 1

    call assert_alloc (work, nrwrk=n*n, niwrk=n)

    shp = n
    call workspace_get_ptr (work, shp, ptr_a)
    call workspace_get_ptr (work, n, ptr_ipiv)

    call vander (x, ptr_a, increasing=.true.)
    ! Copy into COEFS which will contains the coefficients after calling GESV
    coefs(:,:) = y

    ! Fit exactly identified system using GESV
    call GESV (n, nrhs, ptr_a, lda, ptr_ipiv, coefs, ldb, info)
    status%code_orig = info

    ! Map GESV error code to NUMFORT status code
    if (info == 0) then
        status = NF_STATUS_OK
    else if (info < 0) then
        status = NF_STATUS_INVALID_ARG
    else
        status = NF_STATUS_INVALID_STATE
    end if

end subroutine


