
elemental function __APPEND(dgenpareto_pdf,__PREC) (obj, x, loc, scale, shape) &
        result(fx)
    !*  DGENPARETO_PDF implements the PDF for a Generalized Pareto distributed
    !   random variable. 
    !   The PDF is given by 
    !       f(x) = 1/sigma * (1 + xi(x - mu)/sigma) ** (-1/xi - 1)
    !   where mu, sigma and xi are the location, scale and shape parameters,
    !   respectively.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dgenpareto,__PREC)), intent(in) :: obj
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

    lloc = obj%loc
    lscale = obj%scale
    lshape = obj%shape

    if (present(loc)) lloc = loc
    if (present(scale)) lscale = scale
    if (present(shape)) lshape = shape

    fx = 0.0_PREC

    ! Evaluate only if x is in support, ensuring that f(x) is real valued.
    if (lshape >= 0.0) then
        if (x < lloc) return
    else
        if (x > (lloc - lscale/lshape)) return
    end if

    z = (x - lloc) / lscale
    fx = 1.0_PREC/lscale * (1.0_PREC + lshape * z) ** (-1.0_PREC/lshape - 1.0_PREC)

end function


elemental function __APPEND(dgenpareto_cdf,__PREC) (obj, x, loc, scale, shape) &
        result(fx)
    !*  DGENPARETO_CDF implements the CDF for a Generalized Pareto distributed
    !   random variable. 
    !   The CDF is given by 
    !       f(x) = 1 - (1 + xi(x - mu)/sigma) ** (-1/xi)
    !   where mu, sigma and xi are the location, scale and shape parameters,
    !   respectively.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dgenpareto,__PREC)), intent(in) :: obj
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

    lloc = obj%loc
    lscale = obj%scale
    lshape = obj%shape

    if (present(loc)) lloc = loc
    if (present(scale)) lscale = scale
    if (present(shape)) lshape = shape

    fx = 0.0_PREC

    if (lshape >= 0.0) then
        if (x < lloc) return
    else
        ! Support is given by [mu, mu + sigma/xi] according to wiki
        if (x > (lloc - lscale/lshape)) then
            fx = 1.0_PREC
            return
        end if
    end if

    z = (x - lloc) / lscale
    fx = 1.0_PREC - (1.0_PREC + lshape * z) ** (-1.0_PREC/lshape)

end function

