

module numfort_optimize_common

    implicit none

    enum, bind(c)
        enumerator :: OPTIM_PRINT_NONE = 0, &
            OPTIM_PRINT_MINIMAL = 10, &
            OPTIM_PRINT_VERBOSE = 20, &
            OPTIM_PRINT_ALL = 30
    end enum

end module
