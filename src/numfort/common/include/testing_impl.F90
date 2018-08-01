
elemental function __APPEND(is_close,__PREC) (actual, desired, rtol, atol) &
        result(res)
    !*  ALL_CLOSE returns true if two values are close up to a given
    !   tolerance. Specifically, if
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in)  :: actual
    real (PREC), intent(in)  :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    real (PREC) :: lrtol, latol, diff

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    diff = abs(actual-desired)
    res = diff <= (latol + lrtol * abs(desired))

end function


pure function __APPEND(all_close,__PREC) (actual, desired, rtol, atol) &
        result(res)
    !*  ALL_CLOSE returns true if two values are close up to a given
    !   tolerance. Specifically, if
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in)  :: actual
    real (PREC), intent(in)  :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    real (PREC) :: lrtol, latol, diff

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    diff = abs(actual-desired)
    res = diff <= (latol + lrtol * abs(desired))

end function


pure function __APPEND(all_close_1d,__PREC) (actual, desired, rtol, atol) &
        result(res)
    !*  ALL_CLOSE returns true if two values are close up to a given
    !   tolerance. Specifically, if elementwise it holds that
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: actual
    real (PREC), intent(in), dimension(:) :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    real (PREC) :: lrtol, latol, diff
    integer :: i

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    res = .false.
    if (size(actual) /= size(desired)) return

    res = .true.
    do i = 1, size(actual)
        ! Exit immediately if we encounter element pair where
        ! difference between actual and desired exceeds permissible threshold.
        diff = abs(actual(i) - desired(i))
        if (diff > (latol + lrtol * abs(desired(i)))) then
            res = .false.
            return
        end if
    end do

end function


pure function __APPEND(all_close_2d,__PREC) (actual, desired, rtol, atol) result(res)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: actual
    real (PREC), intent(in), dimension(:,:) :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    integer :: m, n
    logical, dimension(:,:), allocatable :: work

    if (.not. shape_equal (actual, desired)) then
        res = .false.
        return
    end if

    m = size(actual, 1)
    n = size(actual, 2)

    allocate (work(m, n))

    work(:,:) = is_close (actual, desired, rtol, atol)
    res = all(work)

    deallocate (work)

end function


pure function __APPEND(all_close_3d,__PREC) (actual, desired, rtol, atol) result(res)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:,:) :: actual
    real (PREC), intent(in), dimension(:,:,:) :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    logical, dimension(:,:,:), allocatable :: work

    if (.not. shape_equal (actual, desired)) then
        res = .false.
        return
    end if

    allocate (work(size(actual,1), size(actual,2), size(actual,3)))

    work(:,:,:) = is_close (actual, desired, rtol, atol)
    res = all(work)

    deallocate (work)

end function
