

subroutine __APPEND(fvec_scalar_args,__PREC) (x, fx, fpx, args)
    import __PREC
    real (__PREC), intent(in), dimension(:) :: x
    real (__PREC), intent(out), optional :: fx
    real (__PREC), intent(out), dimension(:), optional :: fpx
    real (__PREC), intent(in), dimension(:), optional :: args
end subroutine

subroutine __APPEND(fvec_vec_args,__PREC) (x, fx, fpx, args)
    import __PREC
    real (__PREC), intent(in), dimension(:) :: x
    real (__PREC), intent(out), dimension(:), optional :: fx
    real (__PREC), intent(out), dimension(:,:), optional :: fpx
    real (__PREC), intent(in), dimension(:), optional :: args
end subroutine

subroutine __APPEND(fvec_scalar,__PREC) (x, fx, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:) :: x
    real (__PREC), intent(out), optional :: fx
    real (__PREC), intent(out), dimension(:), optional :: fpx
end subroutine

subroutine __APPEND(fvec_vec,__PREC) (x, fx, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:) :: x
    real (__PREC), intent(out), dimension(:), optional :: fx
    real (__PREC), intent(out), dimension(:,:), optional :: fpx
end subroutine
