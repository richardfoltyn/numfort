module numfort

    use iso_fortran_env, only : real64, int64

    implicit none
    private

    public :: linspace, diag, identity, cumsum

    interface linspace
        procedure :: linspace_real64
    end interface

    interface diag
        procedure :: diag_vec_real64, diag_mat_real64
    end interface

    interface identity
        module procedure identity_real64
    end interface

    interface cumsum
        module procedure cumsum_1d_real64, cumsum_2d_real64, cumsum_3d_real64
    end interface

    interface factorial
        module procedure factorial_int
    end interface

    ! Math constants
    real (real64), parameter :: PI = 3.141592653589793238462643383279502884d0

contains

pure function factorial_int(n) result(res)

    integer, intent(in) :: n
    integer :: res
    integer :: i

    res = 1

    do i = 2, n
        res = res * i
    end do

end function

function linspace_real64(from, to, n) result(v)

    real (real64), intent(in) :: from, to
    integer, intent(in) :: n
    real (real64) :: v(n)

    real (real64) :: step
    integer :: i

    step = (to - from) / (n - 1)

    v = [(from + i * step, i=0,n)]
    ! avoid rounding errors in end point
    v(n) = to

end function

subroutine identity_real64(out)

    real (real64), intent(out), dimension(:,:) :: out
    integer :: i

    out = 0.0d0

    forall (i=1:size(out, 1)) out(i,i) = 1.0d0

end subroutine

pure function diag_vec_real64(v) result(res)

    real (real64), intent(in), dimension(:) :: v
    real (real64), dimension(size(v), size(v)) :: res

    integer :: i

    res = 0.0d0

    forall (i=1:size(v)) res(i,i) = v(i)

end function

pure function diag_mat_real64(m) result(res)

    real (real64), intent(in), dimension(:,:) :: m
    real (real64), dimension(size(m, 1)) :: res

    integer :: i

    forall (i=1:size(m, 1)) res(i) = m(i,i)

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

    include "cumsum_wrapper.finc"

end subroutine

subroutine cumsum_3d_real64(x, res, axis)

    real (real64), intent(in), dimension(:,:,:), target, contiguous :: x
    real (real64), dimension(:,:,:), target, contiguous :: res
    real (real64), dimension(:), pointer :: ptr_x => null(), ptr_res => null()

    include "cumsum_wrapper.finc"

end subroutine

subroutine cumsum_nd_real64(x, shp, res, axis)

    real (real64), intent(in), dimension(:) :: x
    real (real64), dimension(:) :: res

    include "cumsum_impl.finc"

end subroutine

subroutine mean_array_2d_real64(data, res, axis)

    real (real64), intent(in), dimension(:,:), target, contiguous :: data
    real (real64), dimension(:), pointer :: ptr_data => null()
    real (real64), intent(out), dimension(:), allocatable :: res

    integer, intent(in), optional :: axis
    integer :: n, laxis, rnk, i
    integer, dimension(:), allocatable :: shp, stride

    if (present(axis)) laxis = axis

    ! rank() does not seem to be in any Fortran standard
    rnk = size(shape(data))

    if (axis > rnk) then
        stop "axis argument exceeds array rank!"
    end if

    allocate (shp(rnk), stride(rnk))
    shp = shape(data)
    n = size(data)

    stride = 1
    do i = 2,rnk
        stride(i) = stride(i-1) * shp(i-1)
    end do

    ptr_data(1:n) => data




    deallocate (shp)

end subroutine


end module
