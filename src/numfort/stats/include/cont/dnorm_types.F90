
type, public :: __APPEND(dnorm,__PREC)
    !*  Type representing Normal distribution.

    real (__PREC) :: loc = 0.0
        !*  Location parameter
    real (__PREC) :: scale = 1.0
        !*  Scale parameter
end type

interface pdf
    module procedure __APPEND(dnorm_pdf,__PREC)
end interface

interface cdf
    module procedure __APPEND(dnorm_cdf,__PREC)
end interface

interface rvs
    module procedure __APPEND(dnorm_rvs,__PREC)
end interface

interface get_dist_params 
    module procedure __APPEND(get_dist_params,__PREC)
end interface
