


module numfort_core_constants
    !*  Module containts common compile-time math constants

    use, intrinsic :: iso_fortran_env

    implicit none

    real (real64), parameter :: PI_real64 = 3.141592653589793238462643383279502884d0
    real (real32), parameter :: PI_real32 = 3.141592653589793238462643383279502884

    real (real64), parameter :: EULER_MASCHERONI_real64 = &
        0.57721566490153286060651209008240243104215933593992_real64
    real (real32), parameter :: EULER_MASCHERONI_real32 = &
        0.57721566490153286060651209008240243104215933593992_real32

    real (real64), parameter :: PI = PI_real64
    real (real64), parameter :: EULER_MASCHERONI = EULER_MASCHERONI_real64


end module
