module numfort_stats_qrng

    use, intrinsic :: iso_fortran_env
    use sobol_direction_num
    use numfort_common

    implicit none

    private

    public :: sobol_init

    integer, parameter :: SOBOL_MAX_NUM_M = 63

    type, public :: sobol_state
        private
        integer :: d = 1
        integer (int64) :: idx = 0
        integer, dimension(SOBOL_MAX_NUM_M) :: v
    end type

contains


pure subroutine sobol_init (self, d, s, a, m_i, status)
    type (sobol_state), intent(in out) :: self
    integer, intent(in) :: d
    integer, intent(in), optional :: s
    integer, intent(in), optional :: a
    integer, intent(in), dimension(:), optional :: m_i
    type (status_t), intent(out), optional :: status

    integer, dimension(SOBOL_MAX_NUM_M) :: lm
    integer :: ls, la, ifrom, ito, k, i, mk
    logical :: user_m

    user_m = present(s) .and. present(a) .and. present(m_i)

    if (user_m) then
        ls = s
        lm(1:ls) = m_i
        la = a
    else
        ls = SOBOL_DEFAULT_S(d)
        la = SOBOL_DEFAULT_A(d)
        ifrom = sum(SOBOL_DEFAULT_S(:d-1)) + 1
        ito = ifrom + ls - 1
        lm(1:ls) = SOBOL_DEFAULT_M(ifrom:ito)
    end if

    do k = ls+1, SOBOL_MAX_NUM_M
        mk = ieor(2 ** ls * lm(k-ls), lm(k-ls))
        do i = 1, ls-1
            mk = ieor(mk, 2 ** i * lm(k-i) * ibtest(la, i))
        end do
        lm(k) = mk
    end do

    do k = 1, SOBOL_MAX_NUM_M
        self%v(k) = lm(k) / (2 ** k)
    end do

end subroutine

pure function ibtest(i, pos) result(res)
    integer, intent(in) :: i
    integer, intent(in) :: pos
    integer :: res

    if (btest(i, pos)) then
        res = 1
    else
        res = 0
    end if

end function

end module
