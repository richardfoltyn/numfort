
elemental function is_close (actual, desired, rtol, atol) result(res)
    !*  ALL_CLOSE returns true if two values are close up to a given
    !   tolerance. Specifically, if
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.

    real (PREC), intent(in)  :: actual
    real (PREC), intent(in)  :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    real (PREC) :: lrtol, latol

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    res = is_close_impl (actual, desired, lrtol, latol)

end function



elemental function is_close_impl (actual, desired, rtol, atol) result(res)
    !*  ALL_CLOSE returns true if two values are close up to a given
    !   tolerance. Specifically, if
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.
    real (PREC), intent(in) :: actual
    real (PREC), intent(in) :: desired
    real (PREC), intent(in) :: rtol, atol
    logical :: res

    real (PREC) :: diff

    diff = abs(actual-desired)
    res = diff <= (atol + rtol * abs(desired))

end function



pure function all_close_0d (actual, desired, rtol, atol) result(res)
    !*  ALL_CLOSE returns true if two values are close up to a given
    !   tolerance. Specifically, if
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.
    real (PREC), intent(in)  :: actual
    real (PREC), intent(in)  :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    real (PREC) :: lrtol, latol

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    res = is_close (actual, desired, lrtol, latol)

end function



pure function all_close_1d (actual, desired, rtol, atol) result(res)
    !*  ALL_CLOSE returns true if two values are close up to a given
    !   tolerance. Specifically, if elementwise it holds that
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.
    real (PREC), intent(in), dimension(:) :: actual
    real (PREC), intent(in), dimension(:) :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    real (PREC) :: lrtol, latol

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    res = .false.
    if (size(actual) /= size(desired)) return

    res = all(is_close (actual, desired, lrtol, latol))

end function



pure function all_close_fast_1d (actual, desired, rtol, atol) result(res)
    !*  ALL_CLOSE_FAST returns true if two values are close up to a given
    !   tolerance. Specifically, if elementwise it holds that
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.
    !
    !   Note: Routine assumes CONTIGUOUS arguments, since some compilers
    !   manage to screw this up despite making a copy at runtime.
    real (PREC), intent(in), dimension(:), contiguous :: actual
    real (PREC), intent(in), dimension(:), contiguous :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    real (PREC) :: lrtol, latol

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    res = .false.
    if (size(actual) /= size(desired)) return

    res = all(is_close (actual, desired, lrtol, latol))

end function



pure function all_close_fast_2d (actual, desired, rtol, atol) result(res)
    !*  ALL_CLOSE_FAST returns true if two values are close up to a given
    !   tolerance. Specifically, if elementwise it holds that
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.
    !
    !   Note: Routine assumes CONTIGUOUS arguments, since some compilers
    !   manage to screw this up despite making a copy at runtime.
    real (PREC), intent(in), dimension(:,:), contiguous :: actual
    real (PREC), intent(in), dimension(:,:), contiguous :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    real (PREC) :: lrtol, latol

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    if (.not. shape_equal (actual, desired)) then
        res = .false.
        return
    end if

    res = all(is_close(actual, desired, lrtol, latol))

end function



pure function all_close_2d (actual, desired, rtol, atol) result(res)
    !*  ALL_CLOSE returns true if two values are close up to a given
    !   tolerance. Specifically, if elementwise it holds that
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.
    real (PREC), intent(in), dimension(:,:) :: actual
    real (PREC), intent(in), dimension(:,:) :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res

    real (PREC) :: lrtol, latol

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    if (.not. shape_equal (actual, desired)) then
        res = .false.
        return
    end if

    res = all(is_close(actual, desired, lrtol, latol))

end function



pure function all_close_3d (actual, desired, rtol, atol) result(res)
    !*  ALL_CLOSE returns true if two values are close up to a given
    !   tolerance. Specifically, if elementwise it holds that
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.
    real (PREC), intent(in), dimension(:,:,:) :: actual
    real (PREC), intent(in), dimension(:,:,:) :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res
    real (PREC) :: lrtol, latol

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    if (.not. shape_equal (actual, desired)) then
        res = .false.
        return
    end if

    res = all(is_close(actual, desired, lrtol, latol))

end function



pure function all_close_fast_3d (actual, desired, rtol, atol) result(res)
    !*  ALL_CLOSE_FAST returns true if two values are close up to a given
    !   tolerance. Specifically, if elementwise it holds that
    !       abs(actual-desired) < atol + rtol * abs(desired)
    !   the function returns true, and false otherwise.
    !
    !   Note: Routine assumes CONTIGUOUS arguments, since some compilers
    !   manage to screw this up despite making a copy at runtime.
    real (PREC), intent(in), dimension(:,:,:), contiguous :: actual
    real (PREC), intent(in), dimension(:,:,:), contiguous :: desired
    real (PREC), intent(in), optional :: rtol, atol
    logical :: res
    real (PREC) :: lrtol, latol

    latol = 0.0_PREC
    lrtol = sqrt(epsilon(0.0_PREC))
    if (present(atol)) latol = atol
    if (present(rtol)) lrtol = rtol

    if (.not. shape_equal (actual, desired)) then
        res = .false.
        return
    end if

    res = all(is_close(actual, desired, lrtol, latol))

end function
