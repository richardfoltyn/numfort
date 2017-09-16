
elemental function __APPEND(dnorm_pdf,__PREC) (obj, x, loc, scale) &
        result(fx)
    !*  DNORM_PDF returns the value of the PDF of a Normal random 
    !   variable at a given point.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dnorm,__PREC)), intent(in) :: obj
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
    real (PREC) :: b

    call get_dist_params (obj, loc, scale, lloc, lscale)

    b = 2 * (lscale ** 2)
    fx = NORM_CONST/lscale * exp(-((x-lloc) ** 2) / b)
end function


impure elemental function __APPEND(dnorm_cdf,__PREC) (obj, x, loc, scale) &
        result(fx)
    !*  DNORM_CDF returns the value of the CDF of a Normal random variable 
    !   at a given point.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dnorm,__PREC)), intent(in) :: obj
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
    ! CDFLIB90 argument
    integer, parameter :: which = 1
    ! CDFLIB90 refuses to operate on domain outside of [-immense, immense]
    real (PREC), parameter :: immense = 1.0e100_PREC

    call get_dist_params (obj, loc, scale, lloc, lscale)

    ! Check whether we are outside of supported domain
    if ((x-lloc) < -immense) then
        fx = 0.0_PREC
    else if ((x-lloc) > immense) then
        fx = 1.0_PREC
    else
        ! call CDFLIB90 routine to actually do the work
        call cdf_normal (which, cum=fx, x=x, mean=lloc, sd=lscale)
    end if

end function


impure elemental subroutine __APPEND(dnorm_rvs,__PREC) (obj, x, loc, scale)
    !*  DNORM_RVS draws a random value from norm distribution with 
    !   given parameters.
    !
    !   Note: Cannot be implemented as ELEMENTAL FUNCTION as then the
    !   the function call is identical for all array elements if 
    !   loc/scale parameters are scalars, and hence the same random value
    !   will be assigned to all elements in result array!

    integer, parameter :: PREC = __PREC

    type (__APPEND(dnorm,__PREC)), intent(in) :: obj
        !*  Object identifying type of distribution
    real (PREC), intent(out) :: x
        !*  Randomly drawn value x.
    real (PREC), intent(in), optional :: loc
        !*  Optional location parameter; if present, overrides
        !   value stored in distribution object
    real (PREC), intent(in), optional :: scale
        !*  Optional scale parameter; if present, overrides
        !   value stored in distribution object

    real (PREC) :: lloc, lscale, z

    call get_dist_params (obj, loc, scale, lloc, lscale)

    ! random draw from std. normal distribution
    z = real(random_normal (), PREC)
    x = lloc + z*lscale
end subroutine


pure subroutine __APPEND(get_dist_params,__PREC) (obj, loc, scale, loc_out, scale_out)
    !*  GET_DIST_PARAMS returns distribution parameters based on
    !   user-provided input and default values.

    integer, parameter :: PREC = __PREC

    type (__APPEND(dnorm,__PREC)), intent(in) :: obj
    real (PREC), intent(in), optional :: loc
    real (PREC), intent(in), optional :: scale
    real (PREC), intent(out) :: loc_out
    real (PREC), intent(out) :: scale_out
 
    loc_out = obj%loc
    scale_out = obj%scale

    if (present(loc)) loc_out = loc
    if (present(scale)) scale_out = scale

end subroutine
