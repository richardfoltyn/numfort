module tests_numfort

    use numfort
    use iso_fortran_env, only : real32, real64, int32, int64

    implicit none

    public :: test_numfort

contains

subroutine test_numfort()

    call test_identity
    call test_diag

end subroutine

subroutine test_identity()

    integer, parameter :: N = 5

    real (real64) :: mat1(N,N), mat2(N,N)

    call identity(mat1)

    where (mat1 == 1.0d0)
        mat2 = 1.0d0
    else where (mat1 == 0.0d0)
        mat2 = 0.0d0
    end where

    if (any(mat1 /= mat2))  then
        stop "Error encountered when creating identity matrix"
    end if

end subroutine


subroutine test_diag()

    integer, parameter :: N = 5

    real (real64), allocatable, dimension(:) :: vec1, vec2
    real (real64), allocatable, dimension(:, :) :: mat1

    vec1 = [1,2,3,4,5]
    mat1 = diag(vec1)
    vec2 = diag(mat1)

    if (any(vec1 /= vec2)) then
        print *, "test_diag: vec1 and vec2 differ"
        print "(a, (f8.4), a)", "vec1=[", vec1, "]"
        print "(a, (f8.4), a)", "vec2=[", vec2, "]"
    end if

end subroutine

end module
