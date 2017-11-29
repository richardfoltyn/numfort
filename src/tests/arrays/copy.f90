program array_copy_test

    use, intrinsic :: iso_fortran_env
    use numfort_arrays
    use numfort_common, only: shape_equal

    use fcore_testing

    implicit none

    integer, parameter :: PREC = real64


    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("numfort_arrays_copy unit tests")

    call test_missing_src (tests)
    call test_unalloc_dst (tests)
    call test_wrong_dst_shape (tests)
    call test_copy (tests)

    call tests%print ()

end subroutine


subroutine test_missing_src (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), allocatable :: src1(:), src2(:,:), src3(:,:,:), src4(:,:,:,:)
    real (PREC), allocatable :: dst1(:), dst2(:,:), dst3(:,:,:), dst4(:,:,:,:)

    tc => tests%add_test ("COPY_ALLOC with missing/unalloc. SRC argument")

    ! 1. Unallocated SRC, unallocated DST argument
    call copy_alloc (src1, dst1)
    call tc%assert_true (.not. allocated(dst1), "1d: unallocated SRC/DST")

    call copy_alloc (src2, dst2)
    call tc%assert_true (.not. allocated(dst2), "2d: unallocated SRC/DST")

    call copy_alloc (src3, dst3)
    call tc%assert_true (.not. allocated(dst3), "3d: unallocated SRC/DST")

    call copy_alloc (src4, dst4)
    call tc%assert_true (.not. allocated(dst4), "4d: unallocated SRC/DST")

    ! 2. Missing SRC argument
    call copy_alloc (dst=dst1)
    call tc%assert_true (.not. allocated(dst1), "1d: missing SRC, unallocated DST")

    call copy_alloc (dst=dst2)
    call tc%assert_true (.not. allocated(dst2), "2d: missing SRC, unallocated DST")

    call copy_alloc (dst=dst3)
    call tc%assert_true (.not. allocated(dst3), "3d: missing SRC, unallocated DST")

    call copy_alloc (dst=dst4)
    call tc%assert_true (.not. allocated(dst4), "4d: missing SRC, unallocated DST")

end subroutine


subroutine test_unalloc_dst (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), allocatable :: src1(:), src2(:,:), src3(:,:,:), src4(:,:,:,:)
    real (PREC), allocatable :: dst1(:), dst2(:,:), dst3(:,:,:), dst4(:,:,:,:)

    tc => tests%add_test ("COPY_ALLOC with unalloc. DST argument")

    allocate (src1(10), source=1.0_PREC)
    call copy_alloc (src1, dst1)
    call tc%assert_true (shape_equal(src1, dst1) .and. all(src1==dst1), &
        "1d: allocated SRC, unalloc. DST")

    allocate (src2(10,10), source=1.0_PREC)
    call copy_alloc (src2, dst2)
    call tc%assert_true (shape_equal(src2, dst2) .and. all(src2==dst2), &
        "2d: allocated SRC, unalloc. DST")

    allocate (src3(10,10,10), source=1.0_PREC)
    call copy_alloc (src3, dst3)
    call tc%assert_true (shape_equal(src3, dst3) .and. all(src3==dst3), &
        "3d: allocated SRC, unalloc. DST")

    allocate (src4(10,10,10,10), source=1.0_PREC)
    call copy_alloc (src4, dst4)
    call tc%assert_true (shape_equal(src4, dst4) .and. all(src4==dst4), &
        "4d: allocated SRC, unalloc. DST")

end subroutine


subroutine test_wrong_dst_shape (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), allocatable :: src1(:), src2(:,:), src3(:,:,:), src4(:,:,:,:)
    real (PREC), allocatable :: dst1(:), dst2(:,:), dst3(:,:,:), dst4(:,:,:,:)

    tc => tests%add_test ("COPY_ALLOC with shape(SRC) != shape(DST) argument")

    allocate (src1(10), source=1.0_PREC)
    allocate (dst1(5))
    call copy_alloc (src1, dst1)
    call tc%assert_true (shape_equal(src1, dst1) .and. all(src1==dst1), &
        "1d: allocated SRC, alloc. DST, shape(SRC) != shape(DST)")

    allocate (src2(10,10), source=1.0_PREC)
    allocate (dst2(10,5))
    call copy_alloc (src2, dst2)
    call tc%assert_true (shape_equal(src2, dst2) .and. all(src2==dst2), &
        "2d: allocated SRC, alloc. DST, shape(SRC) != shape(DST)")

    allocate (src3(10,10,10), source=1.0_PREC)
    allocate (dst3(10,10,1))
    call copy_alloc (src3, dst3)
    call tc%assert_true (shape_equal(src3, dst3) .and. all(src3==dst3), &
        "3d: allocated SRC, alloc. DST, shape(SRC) != shape(DST)")

    allocate (src4(10,10,10,10), source=1.0_PREC)
    allocate (dst4(10,10,10,20))
    call copy_alloc (src4, dst4)
    call tc%assert_true (shape_equal(src4, dst4) .and. all(src4==dst4), &
        "4d: allocated SRC, alloc. DST, shape(SRC) != shape(DST)")

end subroutine


subroutine test_copy (tests)
    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC), allocatable :: src1(:), src2(:,:), src3(:,:,:), src4(:,:,:,:)
    real (PREC), allocatable :: dst1(:), dst2(:,:), dst3(:,:,:), dst4(:,:,:,:)

    tc => tests%add_test ("COPY_ALLOC with conformable SRC/DST arguments")

    allocate (src1(10), source=1.0_PREC)
    allocate (dst1(10), source=0.0_PREC)
    call copy_alloc (src1, dst1)
    call tc%assert_true (shape_equal(src1, dst1) .and. all(src1==dst1), &
        "1d: allocated SRC, allocated DST, shape(SRC) == shape(DST)")

    allocate (src2(10,10), source=1.0_PREC)
    allocate (dst2(10,10), source=0.0_PREC)
    call copy_alloc (src2, dst2)
    call tc%assert_true (shape_equal(src2, dst2) .and. all(src2==dst2), &
        "2d: allocated SRC, allocated DST, shape(SRC) == shape(DST)")

    allocate (src3(10,10,10), source=1.0_PREC)
    allocate (dst3(10,10,10), source=0.0_PREC)
    call copy_alloc (src3, dst3)
    call tc%assert_true (shape_equal(src3, dst3) .and. all(src3==dst3), &
        "3d: allocated SRC, allocated DST, shape(SRC) == shape(DST)")

    allocate (src4(10,10,10,10), source=1.0_PREC)
    allocate (dst4(10,10,10,10), source=0.0_PREC)
    call copy_alloc (src4, dst4)
    call tc%assert_true (shape_equal(src4, dst4) .and. all(src4==dst4), &
        "4d: allocated SRC, allocated DST, shape(SRC) == shape(DST)")

end subroutine


end program
