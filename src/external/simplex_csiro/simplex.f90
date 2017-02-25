! Wrapper module to make available SP and DP implementations under a single
! name.

module simplex_csiro

    use simplex_csiro_real64, only: minim_real64 => minim, func_real64 => functn_if
    use simplex_csiro_real32, only: minim_real32 => minim, func_real32 => functn_if
    implicit none
    private

    interface minim
        procedure minim_real32, minim_real64
    end interface

    public :: minim, func_real32, func_real64

end module
