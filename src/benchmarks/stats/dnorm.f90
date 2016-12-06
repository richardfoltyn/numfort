program benchmark_dnorm

    use iso_fortran_env
    use numfort_arrays
    use numfort_stats, only: norm
    implicit none


    call benchmark_cdf ()

contains

subroutine benchmark_cdf ()

    integer, parameter :: N = 10000
    integer :: i, repeat

    real (real64), dimension(:,:), allocatable :: x, fx
    real (real64), dimension(:), allocatable :: vals

    repeat = 1

    allocate (x(N,N), fx(N,N), vals(N))

    call linspace (vals, -1d5, 1d5)

    do i = 1, N
        x(:, i) = vals
    end do

    do i = 1, repeat
        fx = norm%cdf (x)
    end do


end subroutine

end program
