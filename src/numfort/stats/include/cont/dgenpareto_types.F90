
type, public :: __APPEND(dgenpareto,__PREC)
    !*  Type defines contains for parameters of
    !   Generalized Pareto distribution.

    real (__PREC) :: loc
        !*  Location parameter
    real (__PREC) :: shape
        !*  Shape parameter
    real (__PREC) :: scale
        !*  Scale parameter
end type

interface pdf
    module procedure __APPEND(dgenpareto_pdf,__PREC)
end interface

interface cdf
    module procedure __APPEND(dgenpareto_cdf,__PREC)
end interface
