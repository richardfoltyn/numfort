
type :: __APPEND(dmvnorm,__PREC)
    !*  Type representing multivariate normal distribution.

    real (__PREC), dimension(:), allocatable :: mean
        !*  Vector of means
    real (__PREC), dimension(:), allocatable :: cov
        !*  Variance-Covariance matrix
end type


interface rvs
    procedure __APPEND(dmvnorm_rvs_2d,__PREC)
end interface

interface get_dist_params
    procedure __APPEND(get_dist_params,__PREC)
end interface

interface rvs_check_input
    procedure __APPEND(rvs_check_input,__PREC)
end interface
