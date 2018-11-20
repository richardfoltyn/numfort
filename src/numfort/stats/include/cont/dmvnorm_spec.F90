
type :: __APPEND(dmvnorm,__PREC)
    !*  Type representing multivariate normal distribution.

    real (__PREC), dimension(:), allocatable :: mean
        !*  Vector of means
    real (__PREC), dimension(:,:), allocatable :: cov
        !*  Variance-Covariance matrix
end type


interface rvs
    procedure __APPEND(dmvnorm_rvs_1d,__PREC)
end interface

interface rvs
    procedure __APPEND(dmvnorm_rvs_2d,__PREC)
end interface

interface dist_set_params
    procedure __APPEND(dist_set_params,__PREC)
end interface

interface dist_get_params
    procedure __APPEND(dist_get_params,__PREC)
end interface

interface rvs_check_input
    procedure __APPEND(rvs_check_input,__PREC)
end interface
