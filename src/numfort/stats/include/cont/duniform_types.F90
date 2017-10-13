
type, public :: __APPEND(duniform,__PREC)
    !*  Type representing (continuous) uniform distribution

    real (__PREC) :: loc = 0.0
        !*  Location parameter
    real (__PREC) :: scale = 1.0
        !*  Scale parameter
end type

interface pdf
    module procedure __APPEND(duniform_pdf,__PREC)
end interface

interface cdf
    module procedure __APPEND(duniform_cdf,__PREC)
end interface

interface rvs
    module procedure __APPEND(duniform_rvs,__PREC)
end interface

interface get_dist_params 
    module procedure __APPEND(get_dist_params,__PREC)
end interface
