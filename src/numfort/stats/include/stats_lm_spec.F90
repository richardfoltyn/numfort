

interface ols
    procedure __APPEND(ols_1d,__PREC), __APPEND(ols_2d,__PREC)
end interface

interface ols_check_input
    procedure __APPEND(ols_check_input,__PREC)
end interface

interface ols_get_dims
    procedure __APPEND(ols_get_dims,__PREC)
end interface

interface pca
    procedure __APPEND(pca,__PREC)
end interface

interface pca_check_input
    procedure __APPEND(pca_check_input,__PREC)
end interface

interface pca_get_dims
    procedure __APPEND(pca_get_dims,__PREC)
end interface

interface pcr
    procedure __APPEND(pcr_1d,__PREC), __APPEND(pcr_2d,__PREC)
end interface

interface pcr_check_input
    procedure __APPEND(pcr_check_input,__PREC)
end interface

interface pcr_get_dims
    procedure __APPEND(pcr_get_dims,__PREC)
end interface

interface pcr_pca_get_dims
    procedure __APPEND(pcr_pca_get_dims,__PREC)
end interface

interface pcr_pca_check_input
    procedure __APPEND(pcr_pca_check_input,__PREC)
end interface

interface pcr
    procedure __APPEND(pcr_pca_1d,__PREC), __APPEND(pcr_pca_2d,__PREC)
end interface

interface post_estim
    procedure __APPEND(lm_post_estim,__PREC)
end interface
