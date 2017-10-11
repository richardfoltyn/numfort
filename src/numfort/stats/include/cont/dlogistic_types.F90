
type, public :: __APPEND(dlogistic,__PREC)
    !*  Type representing logistic distribution

    real (__PREC) :: loc
        !*  Location parameter
    real (__PREC) :: scale
        !*  Scale parameter
end type

interface pdf
    module procedure __APPEND(dlogistic_pdf,__PREC)
end interface

interface cdf
    module procedure __APPEND(dlogistic_cdf,__PREC)
end interface

interface ppf
    module procedure __APPEND(dlogistic_ppf,__PREC)
end interface

interface rvs
    module procedure __APPEND(dlogistic_rvs,__PREC)
end interface

interface get_dist_params 
    module procedure __APPEND(get_dist_params,__PREC)
end interface
