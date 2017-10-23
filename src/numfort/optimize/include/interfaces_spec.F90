!   Naming convention:
!   fvs_jac_args defines an interface of a function that
!       1. Maps a vector into a scalar
!       2. Takes an (optional) Jacobian as dummy argument
!       3. Takes (optional) ARGS array as dummy argument

subroutine __APPEND(fvs_jac_args,__PREC) (x, fx, fpx, args)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), optional :: fx
    real (__PREC), intent(out), dimension(:), contiguous, optional :: fpx
    real (__PREC), intent(in), dimension(:), optional :: args
end subroutine

subroutine __APPEND(fvv_jac_args,__PREC) (x, fx, fpx, args)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (__PREC), intent(out), dimension(:,:), contiguous, optional :: fpx
    real (__PREC), intent(in), dimension(:), optional :: args
end subroutine

subroutine __APPEND(fvs_jac,__PREC) (x, fx, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), optional :: fx
    real (__PREC), intent(out), dimension(:), contiguous, optional :: fpx
end subroutine

subroutine __APPEND(fvv_jac,__PREC) (x, fx, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (__PREC), intent(out), dimension(:,:), contiguous, optional :: fpx
end subroutine
