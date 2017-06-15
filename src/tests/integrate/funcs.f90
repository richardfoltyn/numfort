module tests_nf_integrate_funcs

    use, intrinsic :: iso_fortran_env
    implicit none

    integer, private, parameter :: PREC = real64
contains

pure function fcn1 (x) result(fx)
    real (PREC), intent(in) :: x
    real (PREC) :: fx

    fx = x ** (1.0_PREC/4)
end function

pure function fcn2 (x) result(fx)
    real (PREC), intent(in) :: x
    real (PREC) :: fx

    fx = x ** (-2.0_PREC) 
end function

pure function fcn3 (x) result(fx)
    real (PREC), intent(in) :: x
    real (PREC) :: fx

    fx = exp (x)
end function

pure function fcn4 (x) result(fx)
    real (PREC), intent(in) :: x
    real (PREC) :: fx

    fx = max(0.0_PREC, x + 0.05_PREC)
end function

end module
