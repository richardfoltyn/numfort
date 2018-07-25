

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

interface percentile
    procedure  __APPEND(percentile_pmf,__PREC), __APPEND(percentile_pmf_scalar,__PREC)
end interface

interface percentile_pmf_check_input
    procedure __APPEND(percentile_pmf_check_input,__PREC)
end interface
