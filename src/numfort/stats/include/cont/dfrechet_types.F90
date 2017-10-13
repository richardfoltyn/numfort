
type, public :: __APPEND(dfrechet,__PREC)
    !*  Container type to store parameters for Frechet distribution.

    real (__PREC) :: loc = 0.0
        !*  Location parameter
    real (__PREC) :: shape = 0.0
        !*  Shape parameter
    real (__PREC) :: scale = 1.0
        !*  Scale parameter
end type

interface pdf
    module procedure __APPEND(dfrechet_pdf,__PREC)
end interface

interface cdf
    module procedure __APPEND(dfrechet_cdf,__PREC)
end interface

interface ppf
    module procedure __APPEND(dfrechet_ppf,__PREC)
end interface

interface get_dist_params 
    module procedure __APPEND(get_dist_params,__PREC)
end interface
