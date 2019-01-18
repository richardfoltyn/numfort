


elemental function __APPEND(dlognorm_pdf,__PREC) (obj, x, loc, scale) &
        result(fx)
    !*  Log-normal PDF routine
    integer, parameter :: PREC = __PREC
    type (__APPEND(dlognorm,__PREC)), intent(in) :: obj
    real (PREC), intent(in) :: x
        !*  Point at which to evaluate PDF
    real (PREC) :: fx
        !*  Value of PDF at point x
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter of the underlying normal distribution.
        !   If present, overrides value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter of the underlying normal distribution.
        !   If present, overrides value stored in distribution object.

    real (PREC), parameter :: PI = __APPEND(PI,__PREC)
    real (PREC), parameter :: SQRT2_PI = 1.0_PREC / sqrt(2.0_PREC * PI)
    real (PREC) :: lloc, lscale
    real (PREC) :: z

    call get_dist_params (obj, loc, scale, lloc, lscale)

    z = (log(x)- lloc) / lscale
    if (x > 0.0_PREC) then
        fx = SQRT2_PI / lscale / x * exp(-z**2.0_PREC / 2.0_PREC)
    else
        fx = 0.0
    end if

end function



elemental function __APPEND(dlognorm_cdf,__PREC) (obj, x, loc, scale) &
        result(fx)
    !*  Log-normal CDF routine
    integer, parameter :: PREC = __PREC
    type (__APPEND(dlognorm,__PREC)), intent(in) :: obj
    real (PREC), intent(in) :: x
        !*  Point at which to evaluate CDF
    real (PREC) :: fx
        !*  Value of CDF at point x
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter of the underlying normal distribution.
        !   If present, overrides value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter of the underlying normal distribution.
        !   If present, overrides value stored in distribution object.

    real (PREC) :: lloc, lscale
    real (PREC) :: z

    call get_dist_params (obj, loc, scale, lloc, lscale)

    z = (log(x)- lloc) / lscale
    if (x > 0.0_PREC) then
        fx = ndtr (z)
    else
        fx = 0.0_PREC
    end if

end function



elemental function __APPEND(dlognorm_ppf,__PREC) (obj, q, loc, scale) &
        result(x)
    !*  Log-normal PPF (percentile point function; or inverse CDF) returns
    !   the quantile corresponding to a given quantile rank.
    integer, parameter :: PREC = __PREC
    type (__APPEND(dlognorm,__PREC)), intent(in) :: obj
    real (PREC), intent(in) :: q
        !*  Quantile (rank) at which to evaluate inverse CDF
        !   Value must be on the interval [0,1], otherwise NaN is returned.
    real (PREC) :: x
        !*  Value of inverse CDF at Q
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter of the underlying normal distribution.
        !   If present, overrides value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter of the underlying normal distribution.
        !   If present, overrides value stored in distribution object.

    real (PREC) :: lloc, lscale

    call get_dist_params (obj, loc, scale, lloc, lscale)

    if (q >= 0.0_PREC .and. q <= 1.0_PREC) then
        ! Compute quantile for a standard-normal distribution
        x = ndtri (q)
        ! Transform to quantile of underlying normal distribution
        x = lloc + lscale * x
        ! Transform to quantile of log-normal distribution
        x = exp(x)
    else
        x = ieee_value (x, IEEE_QUIET_NAN)
    end if

end function



pure subroutine __APPEND(get_dist_params,__PREC) (obj, loc, scale, loc_out, scale_out)
    !*  GET_DIST_PARAMS returns distribution parameters based on
    !   user-provided input and default values.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dlognorm,__PREC)), intent(in) :: obj
    real (PREC), intent(in), optional :: loc
    real (PREC), intent(in), optional :: scale
    real (PREC), intent(out) :: loc_out
    real (PREC), intent(out) :: scale_out

    loc_out = obj%loc
    scale_out = obj%scale

    if (present(loc)) loc_out = loc
    if (present(scale)) scale_out = scale

end subroutine
