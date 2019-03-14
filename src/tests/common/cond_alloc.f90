

program test_common_cond_alloc
    !*  Module contains unit test routines for the array INSERT routine.

    use, intrinsic :: iso_fortran_env
    use numfort_common, only: shape_equal, cond_alloc

    use fcore_testing
    use fcore_common, only: str, str_array, operator(//)

    implicit none

    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("numfort_common HAS_SHAPE unit tests")

    call test_1d (tests)
    call test_2d (tests)

    call tests%print ()

end subroutine


subroutine test_1d (tests)
    !*  Unit tests for in-place insertion into given argument array.

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: n = 2
    integer, dimension(:), allocatable :: arr, src
    integer :: stat, source
    integer :: i, j
    type (str) :: msg

    tc => tests%add_test ("1d-array tests")

    ! Test with scalar interface
    do j = 0, n
        call cond_alloc (arr, j, stat=stat)
        msg = "Scalar API: unallocated ARR, SHP=" // str_array([j], "i0")
        call tc%assert_true (all(shape(arr)==[j]) .and. stat == 0, msg)

        ! Test with allocated array of correct shape
        call cond_alloc (arr, j, stat=stat)
        msg = "Scalar API: allocated ARR, SHP=" // str_array([j], "i0")
        call tc%assert_true (all(shape(arr)==[j]) .and. stat == -1, msg)

        ! Test with source argument
        deallocate (arr)
        source = 10
        call cond_alloc (arr, j, source=source, stat=stat)
        msg = "Scalar API: unallocated ARR, SHP=" // str_array([j], "i0") &
            // ", SOURCE=" // str(source, "i0")
        call tc%assert_true (all(shape(arr)==[j]) .and. stat == 0 &
            .and. all(arr==source), msg)

        ! Test with source argument, allocated ARR
        call cond_alloc (arr, j, source=source, stat=stat)
        msg = "Scalar API: allocated ARR, SHP=" // str_array([j], "i0") &
            // ", SOURCE=" // str(source, "i0")
        call tc%assert_true (all(shape(arr)==[j]) .and. stat == -1, msg)

        deallocate (arr)
    end do

    ! Test with array interface
    do j = 0, n
        call cond_alloc (arr, [j], source=0, stat=stat)
        msg = "Array API: unallocated ARR, SHP=" // str_array([j], "i0")
        call tc%assert_true (all(shape(arr)==[j]) .and. stat == 0, msg)

        ! Test with allocated array of correct shape
        call cond_alloc (arr, [j], source=0, stat=stat)
        msg = "Array API: allocated ARR, SHP=" // str_array([j], "i0")
        call tc%assert_true (all(shape(arr)==[j]) .and. stat == -1, msg)

        ! Test with array source argument, unallocated ARR
        deallocate (arr)
        allocate (src(j))
        forall (i=1:j) src(i) = i
        call cond_alloc (arr, source=src, stat=stat)
        msg = "Array API: unallocated ARR, SHP=" // str_array([j], "i0") &
            // ", SOURCE=" // str_array(src, "i0")
        call tc%assert_true (all(shape(arr)==[j]) .and. stat == 0 &
            .and. all(arr==src), msg)

        ! Test with array source, allocated ARR
        call cond_alloc (arr, source=src, stat=stat)
        msg = "Array API: allocated ARR, SHP=" // str_array([j], "i0")
        call tc%assert_true (all(shape(arr)==[j]) .and. stat == -1, msg)

        deallocate (arr)
        deallocate (src)

    end do

end subroutine


subroutine test_2d (tests)
    !*  Unit tests for in-place insertion into given argument array.

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, parameter :: n = 2
    integer, dimension(:,:), allocatable :: arr, src
    integer, dimension(2) :: shp
    integer, dimension(:), allocatable :: shp2
    integer :: source, stat
    type (str) :: msg
    integer :: i, j, k, l

    tc => tests%add_test ("2d-array tests")

    do i = 1, n
        do j = 1, n
            shp = [i, j]
            ! Unallocated array
            call cond_alloc (arr, shp, stat=stat)
            msg = "Unallocated ARR, SHP=" // str_array(shp, "i0")
            call tc%assert_true (all(shape(arr)==shp) .and. stat==0, msg)

            ! Allocated array of correct size
            call cond_alloc (arr, shp, stat=stat)
            msg = "Allocated ARR of correct shape, SHP=" // str_array(shp, "i0")
            call tc%assert_true (all(shape(arr)==shp) .and. stat==-1, msg)

            ! Allocated array of wrong shape
            shp = shp + 1
            call cond_alloc (arr, shp, stat=stat)
            msg = "Allocated ARR of incorrect shape, SHP=" // str_array(shp, "i0")
            call tc%assert_true (all(shape(arr)==shp) .and. stat==0, msg)

            ! Unallocated ARR, scalar source argument
            deallocate (arr)
            shp = [i, j]
            source = 10
            call cond_alloc (arr, shp, source=source, stat=stat)
            msg = "Unallocated ARR, SHP=" // str_array(shp, "i0") // ", scalar SOURCE"
            call tc%assert_true (all(shape(arr)==shp) .and. stat==0 .and. &
                all(arr==source), msg)

            ! Allocated ARR of wrong shape, array source argument
            shp = shp + 1
            allocate (src(shp(1),shp(2)), source=0)
            forall (k=1:shp(1),l=1:shp(2)) src(k,l) = 10*l + k
            call cond_alloc (arr, source=src, stat=stat)
            msg = "Allocated ARR or wrong size, SHP=" // str_array(shp, "i0") &
                // ", array SOURCE"
            call tc%assert_true (all(shape(arr)==shp) .and. stat==0 .and. &
                all(arr==src), msg)

            ! Test with non-conformable shape argument
            do k = 0, 3
                if (i == 2) cycle

                allocate (shp2(i), source=0)
                if (allocated(arr)) deallocate (arr)
                call cond_alloc (arr, shp2, stat=stat)

                msg = 'Array API: passing incompatible SHP argument of size ' // str(i)
                call tc%assert_true (stat < 0 .and. .not. allocated(arr), msg)

                deallocate (shp2)
            end do

            if (allocated(arr)) deallocate (arr)

            deallocate (src)
        end do
    end do

end subroutine


end program
