program test_numfort_indexing

    use iso_fortran_env
    use corelib_testing
    use numfort_arrays
    implicit none

    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort.arrays unit tests")

    call test_ind2sub (tests)
    ! call test_sub2ind (tests)

    ! print test statistics
    call tests%print ()

end subroutine


subroutine test_ind2sub (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer (INTSIZE), dimension(:), allocatable :: lin, lin2, shp
    integer (INTSIZE), dimension(:, :), allocatable :: sub, sub2
    integer :: i, n, m

    tc => tests%add_test ("ind2sub test cases")

    ! test conversion of 1-d arrays
    n = 100
    allocate (lin(n), sub(n, 1))
    lin = [(i, i=1,n)]
    sub = -1

    call ind2sub ([n], lin, sub)
    call tc%assert_true (all(lin == sub(:, 1)), &
        "Conversion of 1d linear indices")
    deallocate (lin, sub)

    ! check some simple transformations
    allocate (shp(2))
    shp = [7, 8]
    n = product(shp)
    allocate (lin(n), sub(n, size(shp)), sub2(n, size(shp)))
    lin = [(i, i=1,n)]

    call ind2sub (shp, lin, sub)

    ! manually compute subindices
    sub2(:, 1) = modulo (lin - 1, shp(1)) + 1
    sub2(:, 2) = ceiling (lin / real (shp(1)))

    call tc%assert_true (all(sub(1,:) == [1, 1]), "in2dsub, shp = [7, 8]")
    deallocate (lin, sub, sub2, shp)

    ! check fixed point of ind2sub / sub2ind
    call ind2sub_fp (tc, [11], "sub2ind(ind2sub), shp = [11]")
    call ind2sub_fp (tc, [7, 8], "sub2ind(ind2sub), shp = [7, 8]")
    call ind2sub_fp (tc, [3, 5, 15], "sub2ind(ind2sub), shp = [3, 5, 15]")
    call ind2sub_fp (tc, [3, 5, 15, 7], "sub2ind(ind2sub), shp = [3, 5, 15, 7]")
    call ind2sub_fp (tc, [3, 5, 15, 1, 5], "sub2ind(ind2sub), shp = [3, 5, 15, 1, 5]")
end subroutine

subroutine ind2sub_fp (tc, shp, text)
    class (test_case), intent(in), pointer :: tc
    integer, intent(in), dimension(:) :: shp
    character (len=*), intent(in) :: text

    integer (INTSIZE), dimension(:), allocatable :: lin, lin2
    integer (INTSIZE), dimension(:, :), allocatable :: sub
    integer :: i, n

    n = product(shp)
    allocate (lin(n), lin2(n), sub(n, size(shp)))
    lin = [(i, i=1,n)]

    call ind2sub (shp, lin, sub)
    call sub2ind (shp, sub, lin2)

    call tc%assert_true (all(lin == lin2), text)
end subroutine

end program
