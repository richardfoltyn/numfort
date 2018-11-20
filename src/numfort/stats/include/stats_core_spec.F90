

interface mean
    procedure __APPEND(mean_1d,__PREC), __APPEND(mean_2d,__PREC)
end interface

interface std
    procedure __APPEND(std_1d,__PREC), __APPEND(std_2d,__PREC)
end interface

interface mean_std_check_input
    procedure __APPEND(mean_std_check_input,__PREC)
end interface

interface normalize
    procedure __APPEND(normalize_2d,__PREC)
end interface

interface quantile
    procedure __APPEND(quantile_dispatch,__PREC), __APPEND(quantile_dispatch_scalar,__PREC)
end interface

interface quantile_bins
    procedure  __APPEND(quantile_bins,__PREC)
end interface

interface quantile_discrete
    procedure __APPEND(quantile_discrete,__PREC)
end interface

interface quantile_bins_check_input
    procedure __APPEND(quantile_bins_check_input,__PREC)
end interface

interface quantile_discrete_check_input
    procedure __APPEND(quantile_discrete_check_input,__PREC)
end interface

interface cov
    procedure __APPEND(cov,__PREC)
end interface

interface cov_check_input
    procedure __APPEND(cov_check_input,__PREC)
end interface

interface corrcoef
    procedure __APPEND(corrcoef,__PREC)
end interface
