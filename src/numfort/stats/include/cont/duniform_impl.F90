
elemental function __APPEND(duniform_pdf,__PREC) (obj, x, loc, scale) &
        result(fx)
    !*  DUNIFORM_PDF returns the value of the PDF of a uniform random 
    !   variable at a given point.

    integer, parameter :: PREC = __PREC

    type (__APPEND(duniform,__PREC)), intent(in) :: obj
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

    if (x < lloc .or. x > (lloc + lscale)) then
        fx = 0.0_PREC
    else
        fx = 1.0_PREC / lscale
    end if

end function


elemental function __APPEND(duniform_cdf,__PREC) (obj, x, loc, scale) &
        result(fx)
    !*  DUNIFORM_CDF returns the value of the CDF of a uniform random variable 
    !   at a given point.

    integer, parameter :: PREC = __PREC

    type (__APPEND(duniform,__PREC)), intent(in) :: obj
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

    if (x < lloc) then
        fx = 0.0_PREC
    else if (x > (lloc + lscale)) then
        fx = 1.0_PREC
    else
        fx = (x-1.0_PREC) / lscale
    end if

end function


impure elemental subroutine __APPEND(duniform_rvs,__PREC) (obj, x, loc, scale)
    !*  DUNIFORM_RVS draws a random value from uniform distribution with 
    !   given parameters.
    !
    !   Note: Cannot be implemented as ELEMENTAL FUNCTION as then the
    !   the function call is identical for all array elements if 
    !   loc/scale parameters are scalars, and hence the same random value
    !   will be assigned to all elements in result array!

    integer, parameter :: PREC = __PREC

    type (__APPEND(duniform,__PREC)), intent(in) :: obj
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

    call random_number (x)
    x = x * lscale + lloc
end subroutine


pure subroutine __APPEND(get_dist_params,__PREC) (obj, loc, scale, loc_out, scale_out)
    !*  GET_DIST_PARAMS returns distribution parameters based on
    !   user-provided input and default values.

    integer, parameter :: PREC = __PREC

    type (__APPEND(duniform,__PREC)), intent(in) :: obj
    real (PREC), intent(in), optional :: loc
    real (PREC), intent(in), optional :: scale
    real (PREC), intent(out) :: loc_out
    real (PREC), intent(out) :: scale_out
 
    loc_out = obj%loc
    scale_out = obj%scale

    if (present(loc)) loc_out = loc
    if (present(scale)) scale_out = scale

end subroutine
