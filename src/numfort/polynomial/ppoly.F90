


#include <numfort.h>

module numfort_polynomial_ppoly

    use, intrinsic :: iso_fortran_env

    implicit none
    private

    public :: ppoly_size

    contains

pure function ppoly_size (n, k) result(res)
    !*  PPOLY_SIZE returns the required minimum array size to
    !   store the data of a piecewise polynomial of given degree K with N break
    !   points.
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer :: res

    res = (n-1) * (k + 1)
end function

end module