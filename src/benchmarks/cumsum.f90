program main

    use iso_fortran_env, only : real64
    use numfort

    implicit none

    ! choose dimensions such that CPU cache is too small
    integer, parameter, dimension(3) :: shp = [500, 600, 100]
    integer, parameter :: NREPS = 100

    call benchmark_3d
    ! call benchmark_3d_manual

contains

subroutine fill_3d(x)


    real (real64), dimension(:,:,:) :: x
    integer :: i, j, k, start

    start = 1

    ! populate array block-wise, otherwise expression from array constructor
    ! cannot be allocated in stack
    do i = 1, shp(3)
        do j = 1, shp(2)
            x(:, j, i) = [(k, k=start,start+shp(1))]
            start = start + shp(1)
        end do
    end do

end subroutine

subroutine benchmark_3d()

    real (real64), dimension(:,:,:), allocatable :: x, res

    integer :: i

    allocate (x(shp(1), shp(2), shp(3)))
    allocate (res, mold=x)

    call fill_3d (x)

    do i = 1,NREPS
        call cumsum(x, res, axis=1)
    end do

    deallocate (x, res)

end subroutine

subroutine benchmark_3d_manual()

    real (real64), dimension(:,:,:), allocatable :: x, res

    integer :: i, j, k, l

    allocate (x(shp(1), shp(2), shp(3)))
    allocate (res, mold=x)

    call fill_3d (x)

    do i = 1,NREPS

        res = 0

        do l = 1,shp(3)
            do k = 1,shp(2)
                res(1, k, l) = x(1, k, l)
                do j = 2,shp(1)
                    res(j, k, l) = res(j-1, k, l) + x(j, k, l)
                end do
            end do
        end do
    end do

    deallocate (x, res)

end subroutine

end program main
