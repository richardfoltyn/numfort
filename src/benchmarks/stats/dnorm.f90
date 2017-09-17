program benchmark_dnorm

    use iso_fortran_env
    use numfort_arrays
    use numfort_stats, only: dnorm => dnorm_real64, cdf
    implicit none


    call benchmark_cdf ()

contains

subroutine benchmark_cdf ()

    integer, parameter :: N = 10000
    integer :: i, repeat

    real (real64), dimension(:,:), allocatable :: x, fx
    real (real64), dimension(:), allocatable :: vals
    type (dnorm) :: norm = dnorm(loc=0.0d0, scale=1.0d0)

    repeat = 1

    allocate (x(N,N), fx(N,N), vals(N))

    call linspace (vals, -1d5, 1d5)

    do i = 1, N
        x(:, i) = vals
    end do

    do i = 1, repeat
        fx = cdf (norm, x)
    end do


end subroutine

end program
