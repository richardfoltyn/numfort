! ------------------------------------------------------------------------------
! Input checker for mean and std routines

pure subroutine __APPEND(mean_std_check_input,__PREC) (nvars, nobs, out1, out2, dof, status)
    !*  Check user-supplied input (other than dim) for errors

    integer, parameter :: PREC = __PREC

    integer, intent(in), optional :: nvars
        !!  Number of variables in 2d array
    integer, intent(in) :: nobs
        !!  Number of observations
    real (PREC), intent(in), dimension(:) :: out1
        !!  Array to store mandatory output (mean or std)
    real (PREC), intent(in), dimension(:), optional :: out2
        !!  Optional array to store additional output (mean)
    integer, intent(in), optional :: dof
        !!  Optional degrees of freedom correction
    type (status_t), intent(out) :: status

    status = NF_STATUS_INVALID_ARG

    if (nobs < 1) return

    if (present(nvars)) then
        if (nvars < 1) return
        if (size(out1) < nvars) return
        if (present(out2)) then
            if (size(out2) < nvars) return
        end if
    end if

    if (present(dof)) then
        if (dof < 0) return
    end if

    status = NF_STATUS_OK
end subroutine


! ------------------------------------------------------------------------------
! MEAN for 1d arrays

pure subroutine __APPEND(mean_1d,__PREC) (x, m, dim, status)
    !*  Compute mean of values in 1d array.

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: x
        !!  Array containing numbers for which mean should be computed.
    real (PREC), intent(out) :: m
        !!  Computed mean.
    integer, intent(in), optional :: dim
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays
    type (status_t), intent(out), optional :: status
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays

    ! compute in double precision regardless of input/output kind
    real (PREC) :: m_old, m_new
    type (status_t) :: lstatus
    integer :: i

    ! initialize to default values
    lstatus = NF_STATUS_OK
    m = 0.0_PREC

    if (size(x) == 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! compute running mean using Welford's method; see
    ! http://www.johndcook.com/blog/standard_deviation/

    m_new = x(1)
    m_old = x(1)

    do i = 2, size(x)
        m_new = m_old + (x(i) - m_old) / i
        m_old = m_new
    end do
    m = m_new

100 continue
    if (present(status)) status = lstatus
end subroutine

! ------------------------------------------------------------------------------
! MEAN routine for 2d input arrays

pure subroutine __APPEND(mean_2d,__PREC) (x, m, dim, status)
    !*  MEAN computes the mean of array elements along a dimension.

    integer, parameter :: PREC = __PREC

    ! Input data. Multiple 'variables' are supported
    ! and can be organized either in rows or columns (goverened by the dim argument)
    real (PREC), intent(in), dimension(:,:) :: x
    ! Contains variable means on exit.
    ! Array needs to be at least as large as the number of variables.
    real (PREC), intent(out), dimension(:) :: m
    ! If present, specifies the dimension along which to normalize (default: dim=1)
    integer, intent(in), optional :: dim
    ! If present, returns 0 on exit if normalization was successful,
    ! and status > 0 if an error was encountered
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:), allocatable :: m_old, m_new
    type (status_t) :: lstatus
    ! iv indexes variables, iobs individual observations
    integer :: iv, iobs, nvars, nobs, ldim

    ! default values
    lstatus = NF_STATUS_OK
    m = 0.0_PREC

    call mean_std_init (shape(x), dim, ldim, nvars, nobs, status=lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    call mean_std_check_input (nvars, nobs, m, status=lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! compute mean and std. dev. using Welford's method; see
    ! http://www.johndcook.com/blog/standard_deviation/

    if (ldim == 1) then
        ! If variables are stored in columns, use 1d mean routine
        do iv = 1, nvars
            call mean (x(:, iv), m(iv))
        end do
    else if (ldim == 2) then
        ! Compute running mean and std if variables are stored in rows
        allocate (m_new(nvars), m_old(nvars))
        m_new(:) = x(:, 1)
        m_old(:) = x(:, 1)

        ! Compute running means for all variables simultaneuosly
        do iobs = 2, nobs
            m_new(:) = m_old + (x(:, iobs) - m_old) / iobs
            m_old(:) = m_new
        end do

        m(:) = m_new
    end if

100 continue
    if (present(status)) status = lstatus
end subroutine

! ------------------------------------------------------------------------------
! SD for 1d input arrays

pure subroutine __APPEND(std_1d,__PREC) (x, s, m, dof, dim, status)
    !*  Compute standard deviation (and optionally mean) for data given in
    !   1d array.

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: x
        !   Array containing numbers for which standard deviation should be
        !   computed
    real (PREC), intent(out) :: s
        !   Contains computed standard deviation on exit.
    real (PREC), intent(out), optional :: m
        !   If present, contains computed mean on exit.
    integer, intent(in), optional :: dof
        !*  Degrees of freedom correction to use (default: 1). Computes
        !   standard deviation as sum of squared deviations divided by (N-dof).
    integer, intent(in), optional :: dim
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays
    type (status_t), intent(out), optional :: status
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays

    real (PREC) :: m_old, m_new, s_old, s_new, lm, x_i
    integer :: i, ldof
    type (status_t) :: lstatus

    ! initiliaze default values
    lstatus = NF_STATUS_OK
    ldof = 1
    s = 0.0_PREC
    if (present(m)) m = 0.0_PREC

    if (present(dof)) then
        if (dof < 0) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
        ldof = dof
    end if

    if (size(x) == 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! compute mean and std. dev. using Welford's method; see
    ! http://www.johndcook.com/blog/standard_deviation/

    m_new = x(1)
    m_old = m_new
    s_new = 0.0_PREC
    s_old = 0.0_PREC

    do i = 2, size(x)
        x_i = x(i)
        m_new = m_old + (x_i - m_old) / i
        s_new = s_old + (x_i - m_old) * (x_i - m_new)

        m_old = m_new
        s_old = s_new
    end do

    lm = m_new
    if (size(x) > 1) then
        s = sqrt(s_new / (size(x)-ldof))
    end if

    if (present(m)) m = lm

100 continue
    if (present(status)) status = lstatus
end subroutine

! ------------------------------------------------------------------------------
! SD for 2d input arrays

pure subroutine __APPEND(std_2d,__PREC) (x, s, m, dof, dim, status)
    !*  Compute standard deviation (and optionally) mean of data in 2d-array
    !   along the desired dimension.

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x
        !*  2d-array containing data for which standard deviation should be
        !   computed.

    real (PREC), intent(out), dimension(:) :: s
        !*  On exit, contains standard deviations along desired dimension.
        !   Array size must be sufficient to hold computed values.

    real (PREC), intent(out), dimension(:), optional :: m
        !*  If present, contains vector of means computed along desired dimension.
        !   Array size must be sufficient to hold computed values.
    integer, intent(in), optional :: dof
        !*  Degrees of freedom correction to use (default: 1). Computes
        !   standard deviation as sum of squared deviations divided by (N-dof).
    integer, intent(in), optional :: dim
        !*  If present, specifies dimension along which std/mean should be computed.

    type (status_t), intent(out), optional :: status
        !*  If present, returns NF_STATUS_OK on exit if operation was successful,
        !   and status > 0 if an error was encountered

    real (PREC), dimension(:), allocatable :: lm
    real (PREC), dimension(:), allocatable :: m_old, s_old, s_new
    type (status_t) :: lstatus
    ! iv indexes variables, iobs individual observations
    integer :: iv, iobs, nvars, nobs, ldof, ldim

    ! set default values
    lstatus = NF_STATUS_OK
    s = 0.0_PREC
    ldof = 1
    if (present(m)) m = 0.0_PREC

    call mean_std_init (shape(x), dim, ldim, nvars, nobs, status=lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    call mean_std_check_input (nvars, nobs, m, s, dof, status=lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    if (present(dof)) ldof = dof

    allocate (lm(nvars), source=0.0_PREC)

    ! compute mean and std. dev. using Welford's method; see
    ! http://www.johndcook.com/blog/standard_deviation/

    if (ldim == 1) then
        ! If data is organized in columns we can use 1d routine without
        ! hurting cache locality
        do iv = 1, nvars
            call std (x(:, iv), s(iv), lm(iv), dof=ldof)
        end do
    else if (ldim == 2) then
        ! Compute running mean and std if variables are stored in rows
        ! In this case we don't use 1d method as this would result in
        ! inefficient memory access patterns
        lm(:) = x(:, 1)
        allocate (m_old(nvars), source=x(:, 1))
        allocate (s_new(nvars), source=0.0_PREC)
        allocate (s_old(nvars), source=0.0_PREC)

        ! Compute means / std for all variables simultaneously
        do iobs = 2, nobs
            lm(:) = m_old + (x(:, iobs) - m_old) / iobs
            s_new(:) = s_old + (x(:, iobs) - m_old) * (x(:, iobs) - lm)

            m_old(:) = lm
            s_old(:) = s_new
        end do

        if (nobs > 1) then
            s = sqrt(s_new / (nobs - ldof))
        end if
    end if

    if (present(m)) m = lm

100 continue
    ! Error handling
    if (present(status)) status = lstatus

end subroutine

! ------------------------------------------------------------------------------
! NORMALIZE implementation

pure subroutine __APPEND(normalize_2d,__PREC) (x, m, s, dim, center, scale, &
        dof, status)
    !*  NORMALIZE centers and optionally scales input data
    !   such that the resulting array has zero mean and unity standard deviation
    !   along the specified dimension.

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in out), dimension(:,:) :: x
        !*  Data that should be normalized. Multiple 'variables' are supported
        !   and can be organized either in rows or columns (goverened by
        !   the dim argument)
    real (PREC), intent(out), dimension(:), optional :: m, s
        !*  If present, contain the per-variable means and standard deviations
        !   on exit. Arrays need to be at least as large as the number of variables.
    integer, intent(in), optional :: dim
        !!  If present, specifies the dimension along which to normalize (default: dim=1)
    logical, intent(in), optional :: center
        !!  If present and true, center input data around its mean.
    logical, intent(in), optional :: scale
        !!  If present and true, divide each variable by its standard deviation
    integer, intent(in), optional :: dof
        !*  If present, specifies the degrees-of-freedom correction applied
        !   when computing the standard deviation. Default: 1
    type (status_t), intent(out), optional :: status
        !*  If present, returns 0 on exit if normalization was successful,
        !   and status > 0 if an error was encountered

    real (PREC), dimension(:), allocatable :: lmean, lstd
    integer :: ldim, ldof
    type (status_t) :: lstatus
    ! iv indexes variables, iobs individual observations
    integer :: iv, iobs, nvars, nobs
    logical :: lscale, lcenter

    lcenter = .true.
    lscale = .true.
    ldof = 1
    if (present(center)) lcenter = center
    if (present(scale)) lscale = scale
    if (present(dof)) ldof = dof

    ! initialize dim, nvars and nobs from input data
    call mean_std_init (shape(x), dim, ldim, nvars, nobs, status=lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! Delegate all remaining input checking to STD or MEAN routines

    allocate (lmean(nvars), source=0.0_PREC)
    allocate (lstd(nvars), source=1.0_PREC)

    ! Compute mean (and std) as required.
    ! If only centering is requested (no scaling), compute only
    ! mean and set std = 1.0
    if (lscale) then
        call std (x, s=lstd, m=lmean, dim=ldim, dof=ldof, status=lstatus)
        if (.not. (NF_STATUS_OK .in. lstatus)) goto 100
    else if (lcenter) then
        call mean (x, m=lmean, dim=ldim, status=lstatus)
        if (.not. (NF_STATUS_OK .in. lstatus)) goto 100
    end if

    ! Reset to zero mean correction if case lscale = .true. but lcenter = .false.
    if (.not. lcenter) then
        lmean = 0.0_PREC
    end if

    ! Apply normalization in-place
    if (ldim == 1) then
        ! Compute running mean and std if variables are stored in columns
        do iv = 1, nvars
            ! Apply normalization in-place
            forall (iobs=1:nobs) x(iobs, iv) = (x(iobs, iv) - lmean(iv)) / lstd(iv)
        end do
    else if (ldim == 2) then
        do iobs = 1, nobs
            forall (iv = 1:nvars) x(iv, iobs) = (x(iv, iobs) - lmean(iv)) / lstd(iv)
        end do
    end if

    ! Return mean and std. only if 1) output arrays are present; and 2)
    ! centering / scaling was actually performed.
    if (present(m) .and. lcenter) m = lmean
    if (present(s) .and. lscale) s = lstd

100 continue
    if (present(status)) status = lstatus

end subroutine

! ------------------------------------------------------------------------------
pure subroutine __APPEND(percentile_array,__PREC) (x, q, xq, dim, interp)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: x
        !!  Data array
    real (PREC), intent(in), dimension(:) :: q
        !!  Percentiles to compute which must be between 0 and 100 inclusive.
    real (PREC), intent(out), dimension(:) :: xq
        !!  Output array storing (interpolated) values of x at q
    integer, intent(in), optional :: dim
        !!  Ignored for 1d input arrays, present for compatibility with higher dims.
    character (*), intent(in), optional :: interp
        !!  Interpolation method to use when desired percentile is between
        !!  two data points with indices i and i+1. Valid values are
        !!      1.  'linear': linearly interpolate between values at i, i+1
        !!      2.  'lower': x(i)
        !!      3.  'higher': x(i+1)
        !!      4.  'nearest': value at index which is nearest
        !!      5.  'midpoint': (x(i)+x(i+1))/2

end subroutine

pure subroutine __APPEND(percentile_scalar,__PREC) (x, q, xq, dim, interp)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: x
        !!  Data array
    real (PREC), intent(in) :: q
        !!  Percentile to compute which must be between 0 and 100 inclusive.
    real (PREC), intent(out) :: xq
        !!  (Interpolated) value of x at percentile q
    integer, intent(in), optional :: dim
        !!  Ignored for 1d input arrays, present for compatibility with higher dims.
    character (*), intent(in), optional :: interp
        !!  Interpolation method to use when desired percentile is between
        !!  two data points with indices i and i+1. Valid values are
        !!      1.  'linear': linearly interpolate between values at i, i+1
        !!      2.  'lower': x(i)
        !!      3.  'higher': x(i+1)
        !!      4.  'nearest': value at index which is nearest
        !!      5.  'midpoint': (x(i)+x(i+1))/2

    real (PREC) :: xq1(1), q1(1)

    q1(1) = q
    call percentile (x, q1, xq1, dim, interp)
    xq = xq1(1)

end subroutine
