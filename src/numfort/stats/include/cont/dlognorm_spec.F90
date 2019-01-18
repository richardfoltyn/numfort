
type, public :: __APPEND(dlognorm,__PREC)
    !*  Type representing log-normal distribution where the parameters
    !   LOC and SCALE refer to the underlying Normal distribution.

    real (__PREC) :: loc = 0.0
        !*  Location parameter
    real (__PREC) :: scale = 1.0
        !*  Scale parameter
end type

interface pdf
    procedure __APPEND(dlognorm_pdf,__PREC)
end interface

interface cdf
    procedure __APPEND(dlognorm_cdf,__PREC)
end interface

interface ppf
    procedure __APPEND(dlognorm_ppf,__PREC)
end interface ppf

interface get_dist_params
    procedure __APPEND(get_dist_params,__PREC)
end interface
