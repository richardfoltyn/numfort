
elemental function __APPEND(dlogistic_pdf,__PREC) (obj, x, loc, scale) &
        result(fx)
    !*  DLOGISTIC_PDF returns the value of the PDF of a logistic random 
    !   variable at a given point.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dlogistic,__PREC)), intent(in) :: obj
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

    real (PREC) :: lloc, lscale
    real (PREC) :: ez

    call get_dist_params (obj, loc, scale, lloc, lscale)

    ez = exp( - (x - lloc) / lscale)
    fx = ez / (lscale * (1.0_PREC + ez) ** 2.0_PREC)

end function


elemental function __APPEND(dlogistic_cdf,__PREC) (obj, x, loc, scale) &
        result(fx)
    !*  DLOGISTIC_CDF returns the value of the CDF of a logistic random variable 
    !   at a given point.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dlogistic,__PREC)), intent(in) :: obj
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

    real (PREC) :: lloc, lscale

    call get_dist_params (obj, loc, scale, lloc, lscale)

    fx = 1.0_PREC / (1.0_PREC + exp(- (x-lloc)/lscale))

end function


elemental function __APPEND(dlogistic_ppf,__PREC) (obj, q, loc, scale) &
        result(x)
    !*  DLOGISTIC_PPF (percent point function) returns the value of the 
    !   inverse CDF of a logistic random variable.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dlogistic,__PREC)), intent(in) :: obj
        !*  Object identifying type of distribution
    real (PREC), intent(in) :: q
    real (PREC) :: x
        !*  Value of PDF at point x.
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter; if present, overrides
        !   value stored in distribution object

    real (PREC) :: lloc, lscale

    call get_dist_params (obj, loc, scale, lloc, lscale)

    x = lloc - lscale * log(1.0_PREC / q - 1.0_PREC)
end function


impure elemental subroutine __APPEND(dlogistic_rvs,__PREC) (obj, x, loc, scale)
    !*  DLOGISTIC_RVS draws a random value from logistic distribution with 
    !   given parameters.
    !
    !   Note: Cannot be implemented as ELEMENTAL FUNCTION as then the
    !   the function call is identical for all array elements if 
    !   loc/scale parameters are scalars, and hence the same random value
    !   will be assigned to all elements in result array!

    integer, parameter :: PREC = __PREC

    type (__APPEND(dlogistic,__PREC)), intent(in) :: obj
        !*  Object identifying type of distribution
    real (PREC), intent(out) :: x
        !*  Randomly drawn value x.
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter; if present, overrides
        !   value stored in distribution object

    real (PREC) :: lloc, lscale

    call get_dist_params (obj, loc, scale, lloc, lscale)

    ! TODO: implement PRNG
end subroutine


pure subroutine __APPEND(get_dist_params,__PREC) (obj, loc, scale, loc_out, scale_out)
    !*  GET_DIST_PARAMS returns distribution parameters based on
    !   user-provided input and default values.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dlogistic,__PREC)), intent(in) :: obj
    real (PREC), intent(in), optional :: loc
    real (PREC), intent(in), optional :: scale
    real (PREC), intent(out) :: loc_out
    real (PREC), intent(out) :: scale_out
 
    loc_out = obj%loc
    scale_out = obj%scale

    if (present(loc)) loc_out = loc
    if (present(scale)) scale_out = scale

end subroutine
