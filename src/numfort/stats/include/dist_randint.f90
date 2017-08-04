program test_numfort_stats_randint

    use iso_fortran_env

    use fcore_testing
    use numfort_stats, only: randint, dist_randint

    implicit none

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_stats_randint unit tests")

    call test_rvs (tests)
    ! call test_sub2ind (tests)

    ! print test statistics
    call tests%print ()

end subroutine


subroutine test_rvs (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

#ifdef __SUPPORTS_PDT_KIND__
    type (dist_randint(PREC=real64, INTSIZE=int64)) :: randint64
#endif

    integer (int64), dimension(:), allocatable :: x64
    integer, dimension(:), allocatable :: x
    integer :: i

    tc => tests%add_test ("randint random numbers test cases")

#ifdef __SUPPORTS_PDT_KIND__
    allocate (x64(100))
    call randint64%rvs (x64, low=0_int64, high=10_int64)
#endif

    allocate (x(100))
    call randint%rvs (x)



end subroutine


end program
