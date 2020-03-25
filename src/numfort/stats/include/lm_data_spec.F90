

type :: lm_data
    real (PREC), dimension(:), allocatable :: coefs
        !*  Coefficient array of estimated model
    integer, private :: model = 0
        !*  Type of model estimated (OLS, PCR,...)
    
    real (PREC) :: var_expl
        !*  Frac of RHS variables explained by principle components used 
        !   in PCR (PCR only!)
    logical, private :: add_const = .false.
        !*  True if RHS does not contain intercept and intercept was added 
        !   before estimating the model.
    logical, private :: trans_rhs = .false.
        !*  True if user-provided array of RHS variable needs to be transposed
        !   before estimation.
    integer :: ncomp = 0
        !*  Actual number of principal components used in PCR (PCR only)
    integer :: rank_rhs = 0
        !*  Effective rank of regressor matrix
    integer :: nobs = 0
        !*  Number of observations
    integer :: nvars = 0
        !*  Number of explanatory (RHS) variables 
end type


interface lm_data_update
    procedure lm_data_update
end interface

interface finalize
    procedure lm_data_finalize
end interface

interface assert_alloc_ptr
    procedure lm_data_assert_alloc_ptr
end interface

interface assert_dealloc_ptr
    procedure lm_data_assert_dealloc_ptr
end interface
