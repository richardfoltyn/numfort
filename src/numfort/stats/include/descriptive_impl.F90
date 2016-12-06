! ------------------------------------------------------------------------------
! Input checker for mean and std routines

pure subroutine __APPEND_PREC(mean_std_check_input) (nvars, nobs, out1, out2, dof, status)
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
    integer, intent(out) :: status

    status = STATUS_INVALID_INPUT

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

    status = STATUS_OK
end subroutine


! ------------------------------------------------------------------------------
! MEAN for 1d arrays

pure subroutine __APPEND_PREC(mean_1d) (x, m, dim, status)
    !*  Compute mean of values in 1d array.

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: x
        !!  Array containing numbers for which mean should be computed.
    real (PREC), intent(out) :: m
        !!  Computed mean.
    integer, intent(in), optional :: dim
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays
    integer, intent(out), optional :: status
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays

    ! compute in double precision regardless of input/output kind
    real (real64) :: m_old, m_new
    integer :: i, lstatus

    ! initialize to default values
    lstatus = STATUS_OK
    m = 0.0_PREC

    if (size(x) == 0) then
        lstatus = STATUS_INVALID_INPUT
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

pure subroutine __APPEND_PREC(mean_2d) (x, m, dim, status)
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
    integer, intent(out), optional :: status

    real (real64), dimension(:), allocatable :: m_old, m_new
    integer :: ldim, lstatus
    ! iv indexes variables, iobs individual observations
    integer :: iv, iobs, nvars, nobs

    ! default values
    lstatus = STATUS_OK
    m = 0.0_PREC

    call mean_std_init (shape(x), dim, ldim, nvars, nobs, status=lstatus)
    if (lstatus /= STATUS_OK) goto 100

    call mean_std_check_input (nvars, nobs, m, status=lstatus)
    if (lstatus /= STATUS_OK) goto 100

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

pure subroutine __APPEND_PREC(std_1d) (x, s, m, dof, dim, status)
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
    integer, intent(out), optional :: status
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays

    real (PREC) :: m_old, m_new, s_old, s_new, lm, x_i
    integer :: i, lstatus, ldof

    ! initiliaze default values
    lstatus = STATUS_OK
    ldof = 1
    s = 0.0_PREC
    if (present(m)) m = 0.0_PREC

    if (present(dof)) then
        if (dof < 0) then
            lstatus = STATUS_INVALID_INPUT
            goto 100
        end if
        ldof = dof
    end if
    
    if (size(x) == 0) then
        lstatus = STATUS_INVALID_INPUT
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

pure subroutine __APPEND_PREC(std_2d) (x, s, m, dof, dim, status)
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

    integer, intent(out), optional :: status
        !*  If present, returns STATUS_OK on exit if operation was successful,
        !   and status > 0 if an error was encountered

    real (PREC), dimension(:), allocatable :: lm
    real (PREC), dimension(:), allocatable :: m_old, s_old, s_new
    integer :: ldim, lstatus
    ! iv indexes variables, iobs individual observations
    integer :: iv, iobs, nvars, nobs, ldof

    ! set default values
    status = STATUS_OK
    s = 0.0_PREC
    ldof = 1
    if (present(m)) m = 0.0_PREC

    call mean_std_init (shape(x), dim, ldim, nvars, nobs, status=lstatus)
    if (lstatus /= 0) goto 100

    call mean_std_check_input (nvars, nobs, m, s, dof, status=lstatus)
    if (lstatus /= 0) goto 100

    if (present(dof)) ldof = dof

    allocate (lm(nvars), source=0.0_PREC)

    ! compute mean and std. dev. using Welford's method; see
    ! http://www.johndcook.com/blog/standard_deviation/

    if (ldim == 1) then
        ! If data is organized in columns we can use 1d routine without
        ! hurting cache locality
        do iv = 1, nvars
            call std (x(:, iv), s(iv), lm(iv))
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

100 continue
    if (present(m)) m = lm
    if (present(status)) status = lstatus

end subroutine
