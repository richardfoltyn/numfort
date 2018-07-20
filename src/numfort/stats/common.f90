


module numfort_stats_common

    implicit none
    private

    public :: set_seed

    interface set_seed
        procedure set_seed_1d, set_seed_scalar
    end interface

    contains



subroutine set_seed_1d (val)
    !*  SET_SEED sets the PRNG seed using an array of integer values.
    !   How many elements in VAL are actually used is implementation-dependent,
    !   and depends on the return value N from calling
    !       RANDOM_SEED(size=N)
    !   Hence exactly M=min(size(VAL),N) elements are used to initialize the
    !   seed, and if size(VAL) < N, the remaining elements in the
    !   implementation-specific seed vector are set to 0.
    integer, intent(in), dimension(:) :: val

    integer :: rnd_size, n
    integer, allocatable, dimension(:) :: rnd_put

    call random_seed (size=rnd_size)
    allocate (rnd_put(rnd_size))
    rnd_put(:) = 0

    n = min(size(val), rnd_size)
    rnd_put(1:n) = val(1:n)

    call random_seed (put=rnd_put)

end subroutine


subroutine set_seed_scalar (val)
    !*  SET_SEED sets the PRNG seed using a single scalar value
    integer, intent(in) :: val

    integer, dimension(1) :: val1d

    val1d(1) = val

    call set_seed (val1d)
end subroutine

end module