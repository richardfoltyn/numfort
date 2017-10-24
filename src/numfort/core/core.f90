module numfort_core

    use iso_fortran_env, only : real32, real64, int32, int64


    implicit none
    private

    public :: cumsum, comb, signum
    public :: PI, PI_real32, PI_real64

    interface cumsum
        module procedure cumsum_1d_real64, cumsum_2d_real64, cumsum_3d_real64
    end interface

    interface factorial
        module procedure factorial_int32, factorial_int64
    end interface

    interface comb
        module procedure comb_int32, comb_int64
    end interface

    interface signum
        module procedure signum_real32, signum_real64, signum_int32, signum_int64
    end interface

    ! Math constants
    real (real64), parameter :: PI_real64 = 3.141592653589793238462643383279502884d0
    real (real32), parameter :: PI_real32 = 3.141592653589793238462643383279502884

    real (real64), parameter :: PI = PI_real64

contains

elemental function factorial_int32(n) result(res)
    integer, parameter :: PREC = int32
    integer (PREC), intent(in) :: n
    integer (PREC) :: res
    integer (PREC) :: i

    res = 1
    do i = 2, n
        res = res * i
    end do
end function

elemental function factorial_int64(n) result(res)
    integer, parameter :: PREC = int64
    integer (PREC), intent(in) :: n
    integer (PREC) :: res
    integer (PREC) :: i

    res = 1
    do i = 2, n
        res = res * i
    end do
end function

elemental function comb_int32 (n, k, repetition) result(res)
    integer, parameter :: PREC = int32
    integer (PREC), intent(in) :: n, k
    logical, intent(in), optional :: repetition
    integer (PREC) :: res

    res = int(comb(int(n, int64), int(k, int64), repetition), int32)
end function

elemental function comb_int64 (n, k, repetition) result(res)
    integer, parameter :: PREC = int64
    integer (PREC), intent(in) :: n, k
    logical, intent(in), optional :: repetition
    integer (PREC) :: res

    integer (PREC) :: i

    logical :: lrep
    lrep = .false.

    if (present(repetition)) lrep = repetition

    if (lrep) then
        ! combination with repetition:
        ! this is (n + k - 1)!/(k! (n-1)!)
        ! First compute (n + k - 1)!/(n-1)! = (n+k-1) * ... * n
        res = n
        do i = res+1, n+k-1
            res = res * i
        end do
    else
        ! combination without repetition:
        ! this is n-choose-k, ie \binomial{n}{k}
        ! compute n! / (n-k)! = n * (n-1) * ... * (n-k+1)
        res = n - k + 1
        do i = res + 1, n
            res = res * i
        end do
    end if

    ! Adjust for k! ways to order k elements
    res = res / factorial (k)
end function


!-------------------------------------------------------------------------------
! SIGNUM: alternative SIGN function
! Implement behavior that corresponds to what is found in numpy.

elemental function signum_real64 (x) result(res)
    !*  SIGNUM returns -1, 0 or 1 depending on whether x < 0, x == 0 or x > 0.
    !   Note: Differs from intrinsic SIGN, as SIGN return 1 whenever x >= 0
    !   and -1 otherwise.
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: x
    real (PREC) :: res

    if (x > 0.0_PREC) then
        res = 1.0_PREC
    else if (x < 0.0_PREC) then
        res = -1.0_PREC
    else
        res = 0.0_PREC
    end if
end function


elemental function signum_real32 (x) result(res)
    !*  SIGNUM returns -1, 0 or 1 depending on whether x < 0, x == 0 or x > 0.
    !   Note: Differs from intrinsic SIGN, as SIGN return 1 whenever x >= 0
    !   and -1 otherwise.
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: x
    real (PREC) :: res

    if (x > 0.0_PREC) then
        res = 1.0_PREC
    else if (x < 0.0_PREC) then
        res = -1.0_PREC
    else
        res = 0.0_PREC
    end if
end function


elemental function signum_int32 (x) result(res)
    !*  SIGNUM returns -1, 0 or 1 depending on whether x < 0, x == 0 or x > 0.
    !   Note: Differs from intrinsic SIGN, as SIGN return 1 whenever x >= 0
    !   and -1 otherwise.
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: x
    integer (INTSIZE)  :: res

    if (x > 0_INTSIZE) then
        res = 1_INTSIZE
    else if (x < 0_INTSIZE) then
        res = -1_INTSIZE
    else
        res = 0_INTSIZE
    end if
end function


elemental function signum_int64 (x) result(res)
    !*  SIGNUM returns -1, 0 or 1 depending on whether x < 0, x == 0 or x > 0.
    !   Note: Differs from intrinsic SIGN, as SIGN return 1 whenever x >= 0
    !   and -1 otherwise.
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in) :: x
    integer (INTSIZE)  :: res

    if (x > 0_INTSIZE) then
        res = 1_INTSIZE
    else if (x < 0_INTSIZE) then
        res = -1_INTSIZE
    else
        res = 0_INTSIZE
    end if
end function


!-------------------------------------------------------------------------------
subroutine cumsum_1d_real64(x, res, axis)

    real (real64), intent(in), dimension(:) :: x
    real (real64), dimension(:) :: res
    ! ignore axis for 1D, only relevant for higher-dimensional arrays
    integer, intent(in), optional :: axis

    integer :: i

    res = x(1)

    do i = 2, size(x)
        res(i) = res(i-1) + x(i)
    end do

end subroutine


subroutine cumsum_2d_real64(x, res, axis)

    real (real64), intent(in), dimension(:,:), target, contiguous :: x
    real (real64), dimension(:, :), target, contiguous :: res
    real (real64), dimension(:), pointer :: ptr_x => null(), ptr_res => null()

    include "include/cumsum_wrapper.f90"

end subroutine

subroutine cumsum_3d_real64(x, res, axis)

    real (real64), intent(in), dimension(:,:,:), target, contiguous :: x
    real (real64), dimension(:,:,:), target, contiguous :: res
    real (real64), dimension(:), pointer :: ptr_x => null(), ptr_res => null()

    include "include/cumsum_wrapper.f90"

end subroutine

subroutine cumsum_nd_real64(x, shp, res, axis)

    real (real64), intent(in), dimension(:) :: x
    real (real64), dimension(:) :: res

    include "include/cumsum_impl.f90"

end subroutine


end module
