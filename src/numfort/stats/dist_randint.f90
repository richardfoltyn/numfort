module numfort_stats_dist_randint

    use numfort_stats_disc_dist

    implicit none
    private


    type, extends(disc_dist) :: dist_randint
        integer (int32) :: low
        integer (int32) :: high
    contains
        procedure, pass
    end type

end module
