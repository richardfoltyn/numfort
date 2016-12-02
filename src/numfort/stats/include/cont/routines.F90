! ------------------------------------------------------------------------------
! PDF Method implementation

impure elemental subroutine __APPEND_PREC(pdf) (self, x, fx)
    integer, parameter :: PREC = __PREC
    class (cont_dist __PDT_PARAM_DECL(PREC)), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%pdf (x1, fx1)
    fx = fx1(1)
end subroutine

! ------------------------------------------------------------------------------
! CDF method implementation

subroutine __APPEND_PREC(cdf_scalar) (self, x, fx)
    integer, parameter :: PREC = __PREC
    class (cont_dist __PDT_PARAM_DECL(PREC)), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx

    real (PREC), dimension(1) :: x1, fx1
    x1(1) = x
    call self%cdf (x1, fx1)
    fx = fx1(1)
end subroutine

! ------------------------------------------------------------------------------
! RVS method implementation

subroutine __APPEND_PREC(rvs_scalar) (self, x)
    integer, parameter :: PREC = __PREC
    class (cont_dist __PDT_PARAM_DECL(PREC)), intent(in) :: self
    real (PREC), intent(out) :: x

    real (PREC), dimension(1) :: x1
    x1(1) = x
    call self%rvs (x1)
    x = x1(1)
end subroutine
