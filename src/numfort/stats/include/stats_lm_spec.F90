

integer, parameter :: DEFAULT_PCR_NCOMP = -1

type :: lm_config
    logical :: trans_x = .false.
        !*  Tranpose regressor matrix
    logical :: trans_y = .false.
        !*  Transpose matrix of outcome (LHS) variables. Only relevant
        !   for estimating mutliple outcome variables at once.
    logical :: trans_coefs = .false.
        !*  Transpose coefficient array
    logical :: add_intercept = .true.
        !*  If true, prepend intercept(s) as the first element (or row)
        !   in the coefficient array.
    logical :: scale = .true.
    logical :: center = .true.
    real (PREC) :: var_rhs_min = 0.0_PREC
        !*  Min. fraction of RHS variance captures by principal components
        !   used in PCR.
    integer :: ncomp_min = 0
        !*  Min. number of princal components to use for PCR
    integer :: ncomp_max = -1
    integer :: ncomp = DEFAULT_PCR_NCOMP
        !*  Number of principal components to use for PCR. Overrides any
        !   other setting.
    real (PREC) :: rcond = 100.0_PREC * epsilon(0.0_PREC)
        !*  Argument passed to GELSD for OLS that allows to control the
        !   effective rank of the regressor matrix.
        !   Use same default as Scipy.

    integer :: cv_n = 20
    logical :: cv_parsimonious = .false.
    real (PREC) :: cv_se_mult = 1.0
    logical :: cv_shuffle_obs = .false.
end type

type :: lm_result

    type (lm_config) :: conf

    real (PREC), dimension(:), allocatable :: coefs
        !*  Coefficient array of estimated model
    real (PREC) :: intercept = 0.0_PREC
    real (PREC), dimension(:,:), allocatable :: coefs_multi
    real (PREC), dimension(:), allocatable :: intercept_multi

    integer, private :: model = 0
        !*  Type of model estimated (OLS, PCR,...)

    real (PREC) :: var_rhs = 1.0
        !*  Fraction of RHS variables included in model, e.g. when running
        !   principal component regression.
    integer :: ncomp = 0
        !*  Actual number of principal components used in PCR (PCR only)
    integer :: rank_rhs = 0
        !*  Effective rank of regressor matrix
    integer :: nobs = 0
        !*  Number of observations
    integer :: nrhs = 0
        !*  Number of explanatory (RHS) variables
    integer :: nlhs = 0
        !*  Number of response (dependent, LHS) variables
end type


interface lm_result_update
    procedure lm_result_update
end interface

interface lm_result_reset
    procedure lm_result_reset
end interface

interface ols
    procedure ols_1d, ols_2d
end interface

interface pca
    procedure pca
end interface

interface pcr
    procedure pcr_1d, pcr_2d
end interface

interface pcr
    procedure pcr_pca_1d, pcr_pca_2d
end interface

interface pcr
    procedure pcr_pca_masked_2d
end interface

interface pcr_cv
    procedure pcr_cv_1d, pcr_cv_2d
end interface

interface post_estim
    procedure lm_post_estim
end interface

interface get_dims
    procedure get_dims_2d
end interface

interface predict
    procedure predict_1d, predict_2d, lm_predict_1d, lm_predict_2d
end interface
