
pure subroutine __APPEND(sobol_next_array,__PREC) (self, x, status)
    integer, parameter :: PREC  = __PREC

    type (sobol_state), intent(in out) :: self
    real (PREC), intent(out), dimension(:) :: x
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer (int64) :: k, i

    lstatus = NF_STATUS_OK

    if (self%ndim /= size(x)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (self%idx == 0) then
        ! First point in sequence
        allocate (self%x(self%ndim))
        self%x = 0
        x = 0.0_PREC
        ! Sequence index is 1-based. Note that in the brief documentation
        ! by Joe/Kuo this index is 0-based.
        self%idx = 1
    else
        ! Find first 0 digit from the right in binary representation of
        ! sequence index
        i = self%idx - 1
        do k = 1, SOBOL_MAX_BITS
            if (.not. btest(i, 0)) exit
            i = ishft(i, - 1_int64)
        end do

        self%x(:) = ieor(self%x, self%v(:,k))
        ! Actual x
        x = real(self%x, real64) / (2.0_real64 ** SOBOL_MAX_BITS)

        self%idx = self%idx + 1
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine

pure subroutine __APPEND(sobol_next_scalar,__PREC) (self, x, status)
    integer, parameter :: PREC  = __PREC
    type (sobol_state), intent(in out) :: self
    real (PREC), intent(out) :: x
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(1) :: x1

    x1 = 0.0

    call sobol_next (self, x1, status)

    x = x1(1)

end subroutine
