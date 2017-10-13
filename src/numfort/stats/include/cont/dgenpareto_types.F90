
type, public :: __APPEND(dgenpareto,__PREC)
    !*  Container type to store parameters for Generalized Pareto distribution.

    real (__PREC) :: loc = 0.0
        !*  Location parameter
    real (__PREC) :: shape
        !*  Shape parameter
    real (__PREC) :: scale = 1.0
        !*  Scale parameter
end type

interface pdf
    module procedure __APPEND(dgenpareto_pdf,__PREC)
end interface

interface cdf
    module procedure __APPEND(dgenpareto_cdf,__PREC)
end interface
