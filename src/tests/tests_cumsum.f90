module tests_cumsum

    use iso_fortran_env
    use numfort, only : cumsum

    implicit none
    private

    public :: test_cumsum

contains

subroutine test_cumsum()

    call test_cumsum_1d
    call test_cumsum_2d
    call test_cumsum_3d

end subroutine

subroutine test_cumsum_1d()

    real (real64), dimension(:), allocatable :: x, res, resm

    integer, parameter :: shapes(3) = [1, 10, 100]
    integer :: im, j, m

    do im = 1, size(shapes)
        m = shapes(im)

        allocate (x(m), res(m), resm(m))
        x = [(j, j=1,size(x))]

        call cumsum(x, res)

        resm(1) = x(1)
        do j = 2, size(x)
            resm(j) = resm(j-1) + x(j)
        end do

        if (any(abs(res - resm) > 1d-12)) then
            print *, "CUMSUM[1D] failed; M =", m
            stop
        end if

        deallocate (x, res, resm)

    end do

end subroutine

subroutine test_cumsum_2d()

    real (real64), dimension(:, :), allocatable :: x, res, resm

    integer, parameter :: shapes(3) = [1, 10, 100]
    integer :: i, n, m, im, in

    do im = 1, size(shapes)
        do in = 1, size(shapes)

            m = shapes(im)
            n = shapes(in)

            allocate (x(m, n), res(m, n), resm(m, n))

            x = reshape([(i, i = 1, n*m)], shape=[m, n])

            ! test axis 1
            call cumsum(x, res, axis=1)

            resm(1, :) = x(1, :)
            do i = 2, m
                resm(i, :) = resm(i - 1, :) + x(i, :)
            end do

            if (any(abs(res - resm) > 1d-12)) then
                print *, "CUMSUM[2D] failed; axis=1 M=", m, "N=", n
            end if

            ! text axis 2
            call cumsum(x, res, axis=2)

            resm(:, 1) = x(:, 1)
            do i = 2, n
                resm(:, i) = resm(:, i - 1) + x(:, i)
            end do

            if (any(abs(res - resm) > 1d-12)) then
                print *, "CUMSUM[2D] failed; axis=2 M=", m, "N=", n
            end if

            deallocate (x, res, resm)

        end do
    end do

end subroutine

subroutine test_cumsum_3d()

    real (real64), dimension(:,:,:), allocatable :: x, res, resm

    integer :: im, in, ik, i, j, axis
    integer, parameter :: axis_len(3) = [1, 10, 100]
    integer, dimension(:, :), allocatable :: shapes
    integer, dimension(3) :: shp

    allocate (shapes(size(axis_len), size(axis_len) ** 3))

    ! construct grid of all possible shapes formen from axis lengths specified
    ! above.
    i = 1
    do im = 1, size(axis_len)
        do in = 1, size(axis_len)
            do ik = 1, size(axis_len)
                shapes(:, i) = [axis_len(im), axis_len(in), axis_len(ik)]
                i = i + 1
            end do
        end do
    end do

    do i = 1, size(shapes, 2)
        shp = shapes(:, i)

        allocate (x(shp(1), shp(2), shp(3)))
        allocate (res, mold=x)
        allocate (resm, mold=x)

        ! might exceed stack size
        x = reshape([(j, j=1,product(shp))], shape=shp)

        ! test axis=1
        axis = 1
        call cumsum(x, res, axis=axis)

        resm = 0
        do j = 1, shp(axis)
            resm(j, :, :) = sum(x(1:j, :, :), dim=axis)
        end do

        if (any(abs(res - resm) > 1e-12)) then
            print *, "CUMSUM[3D] failed for axis=", axis, "shape=", shp
        end if

        ! text axis = 2
        axis = 2
        call cumsum(x, res, axis=axis)

        resm = 0
        do j = 1, shp(axis)
            resm(:, j, :) = sum(x(:, 1:j, :), dim=axis)
        end do

        if (any(abs(res - resm) > 1e-12)) then
            print *, "CUMSUM[3D] failed for axis=", axis, "shape=", shp
        end if

        ! text axis = 3
        axis = 3
        call cumsum(x, res, axis=axis)

        resm = 0
        do j = 1, shp(axis)
            resm(:, :, j) = sum(x(:, :, 1:j), dim=axis)
        end do

        if (any(abs(res - resm) > 1e-12)) then
            print *, "CUMSUM[3D] failed for axis=", axis, "shape=", shp
        end if

        deallocate (x)
        deallocate (res)
        deallocate (resm)

    end do

end subroutine


end
