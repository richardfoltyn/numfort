program test_numfort_stats_randint

    use iso_fortran_env

    use fcore_testing
    use numfort_stats, only: randint, drandint

    implicit none

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_stats_randint unit tests")

    call test_rvs (tests)

    ! print test statistics
    call tests%print ()

end subroutine


subroutine test_rvs (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    integer, dimension(:), allocatable :: x
    integer :: i

    tc => tests%add_test ("randint random numbers test cases")

    allocate (x(100))
    call randint%rvs (x)
    
end subroutine


end program
