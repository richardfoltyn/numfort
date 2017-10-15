
elemental function __APPEND(dfrechet_mean,__PREC) (obj, loc, scale, shape) &
        result(fx)
    !*  DFRECHET_MEAN returns the unconditional mean of a Frechet-distributed
    !   r.v.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dfrechet,__PREC)), intent(in) :: obj
        !*  Object identifying type of distribution
    real (PREC) :: fx
        !*  Unconditional mean
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: shape
        !*  Optional shape parameter; if present, overrides
        !   value stored in distribution object

    real (PREC) :: lloc, lscale, lshape

    call get_dist_params (obj, loc, scale, shape, lloc, lscale, lshape)

    if (lshape > 1.0_PREC) then
        fx = lloc + lscale * gamma(1.0_PREC - 1.0_PREC/lshape)
    else
        fx = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    end if
end function


elemental function __APPEND(dfrechet_pdf,__PREC) (obj, x, loc, scale, shape) &
        result(fx)
    !*  DFRECHET_PDF implements the PDF for a Frechet-distributed
    !   random variable.
    !   The PDF is given by
    !       f(x) = alpha/s * ((x-m)/s)**(-1-alpha) * exp(-((x-m)/s)**(-alpha))
    !   where m, s and alpha are the location, scale and shape parameters,
    !   respectively.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dfrechet,__PREC)), intent(in) :: obj
        !*  Object identifying type of distribution
    real (PREC), intent(in) :: x
    real (PREC) :: fx
        !*  Value of PDF at point x.
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: shape
        !*  Optional shape parameter; if present, overrides
        !   value stored in distribution object

    real (PREC) :: lloc, lscale, lshape, z

    call get_dist_params (obj, loc, scale, shape, lloc, lscale, lshape)

    fx = 0.0_PREC

    ! Evaluate only if x is in support, ensuring that f(x) is real valued.
    if (x > lloc) then
        z = (x - lloc) / lscale
        fx = lshape/lscale * z ** (-1.0_PREC-lshape) * exp(-z**(-lshape))
    end if

end function


elemental function __APPEND(dfrechet_cdf,__PREC) (obj, x, loc, scale, shape) &
        result(fx)
    !*  DFRECHET_CDF implements the CDF for a Frechet-distributed
    !   random variable.
    !   The CDF is given by
    !       f(x) = exp(-((x-m)/s) ** (-alpha))
    !   where m, s and alpha are the location, scale and shape parameters,
    !   respectively.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dfrechet,__PREC)), intent(in) :: obj
        !*  Object identifying type of distribution
    real (PREC), intent(in) :: x
    real (PREC) :: fx
        !*  Value of PDF at point x.
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: shape
        !*  Optional shape parameter; if present, overrides
        !   value stored in distribution object

    real (PREC) :: lloc, lscale, lshape, z

    call get_dist_params (obj, loc, scale, shape, lloc, lscale, lshape)

    fx = 0.0_PREC

    if (x > lloc) then
        z = (x - lloc) / lscale
        fx =  exp(- z**(-lshape))
    end if

end function


elemental function __APPEND(dfrechet_ppf,__PREC) (obj, q, loc, scale, shape) &
        result(x)
    !*  DFRECHET_PPF implements the percent point function (inverse CDF)
    !   for a Frechet-distributed random variable.
    !   The PPF is given by
    !       F^{-1}(q) = exp(-((q-m)/s) ** (-alpha))
    !   where m, s and alpha are the location, scale and shape parameters,
    !   respectively.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dfrechet,__PREC)), intent(in) :: obj
        !*  Object identifying type of distribution
    real (PREC), intent(in) :: q
        !*  Percentile at which to compute PPF
    real (PREC) :: x
        !*  Value of PPF at q.
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: shape
        !*  Optional shape parameter; if present, overrides
        !   value stored in distribution object

    real (PREC) :: lloc, lscale, lshape

    call get_dist_params (obj, loc, scale, shape, lloc, lscale, lshape)

    x = lloc

    if (q > 0.0_PREC .and. q < 1.0_PREC) then
        x = lloc + lscale * (-log(q)) ** (-1.0_PREC/lshape)
    else if (q == 1.0_PREC) then
        x = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    end if

end function

pure subroutine __APPEND(get_dist_params,__PREC) (obj, loc, scl, shp, &
        loc_out, scl_out, shp_out)
    !*  GET_DIST_PARAMS returns distribution parameters based on
    !   user-provided input and default values.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dfrechet,__PREC)), intent(in) :: obj
    real (PREC), intent(in), optional :: loc
    real (PREC), intent(in), optional :: scl
    real (PREC), intent(in), optional :: shp
    real (PREC), intent(out) :: loc_out
    real (PREC), intent(out) :: scl_out
    real (PREC), intent(out) :: shp_out

    loc_out = obj%loc
    scl_out = obj%scale
    shp_out = obj%shape

    if (present(loc)) loc_out = loc
    if (present(scl)) scl_out = scl
    if (present(shp)) shp_out = shp

end subroutine
