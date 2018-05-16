

pure subroutine __APPEND(polyfit_check_input,__PREC) (x, y, deg, coefs, status)
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

    real (PREC), dimension(:,:), pointer :: ptr_y
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs

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
    integer, dimension(:), pointer, contiguous :: ptr_ipiv

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


pure subroutine __APPEND(polyfit_deriv_check_input,__PREC) (y, coefs, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(in), dimension(:) :: coefs
    type (status_t), intent(out) :: status

    status = NF_STATUS_INVALID_ARG
    
    if (size(coefs) < size(y)) return
    
    status = NF_STATUS_OK
end subroutine


pure subroutine __APPEND(polyfit_deriv_scalar,__PREC) (x, y, coefs, work, status)
    !*  POLYFIT_DERIV_SCALAR fits a polynomial to a set of function values
    !   and derivatives evaluated at the same given point X.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
        !*  Point X where function value and derivatives are evaluated
    real (PREC), intent(in), dimension(0:) :: y
        !*  List of function values ("order 0" derivative) and k-th order
        !   derivatives in increasing order, ie
        !   y = [f(x), df(x)/x, d^2f(x)/dx^2, ... , d^kf(x)/dx^k]
    real (PREC), intent(out), dimension(0:) :: coefs
        !*  Coefficients of fitted polynomial
    type (__APPEND(workspace,__PREC)), intent(inout), optional :: work
        !*  Optional workspace object
    type (status_t), intent(out), optional :: status
        !*  Optional exit code
    
    type (__APPEND(workspace,__PREC)), pointer :: ptr_work
    real (PREC), dimension(:), pointer, contiguous :: xp, fact
    type (status_t) :: lstatus
    integer :: deg, n, i, j, k, nrwrk
    real (PREC) :: z
    
    lstatus = NF_STATUS_OK
    call polyfit_deriv_check_input (y, coefs, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100
    
    call assert_alloc_ptr (work, ptr_work)
    
    deg = size(y) - 1
    n = deg
    nrwrk = 2 * (n + 1)
    
    call assert_alloc (ptr_work, nrwrk=nrwrk)
    
    fact(0:n) => ptr_work%rwrk(1:n+1)
    xp(0:n) => ptr_work%rwrk(n+2:2*(n+1))
    
    fact(0) = 1.0_PREC
    xp(0) = 1.0_PREC
    do i = 1, n
        fact(i) = fact(i-1) * i
        xp(i) = xp(i-1) * x
    end do
    
    do i = n, 0, -1
        z = y(i)
        do j = n, i+1, -1
            k = j - i
            z = z - fact(j) / fact(k) * coefs(j) * xp(k)
        end do
        
        coefs(i) = z / fact(i)
    end do
    
    
100 continue

    call assert_dealloc_ptr (work, ptr_work)
    if (present(status)) status = lstatus

end subroutine

