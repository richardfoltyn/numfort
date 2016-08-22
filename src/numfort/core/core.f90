module numfort_core

    use iso_fortran_env, only : real32, real64, int32, int64


    implicit none
    private

    public :: cumsum, comb
    public :: PI, PI_REAL32

    interface cumsum
        module procedure cumsum_1d_real64, cumsum_2d_real64, cumsum_3d_real64
    end interface

    interface factorial
        module procedure factorial_int32, factorial_int64
    end interface

    interface comb
        module procedure comb_int32, comb_int64
    end interface

    ! Math constants
    real (real64), parameter :: PI = 3.141592653589793238462643383279502884d0
    real (real32), parameter :: PI_REAL32 = 3.141592653589793238462643383279502884

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
