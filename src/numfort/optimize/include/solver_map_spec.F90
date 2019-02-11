

interface solver_map_init
    procedure __APPEND(solver_map_init,__PREC)
end interface

interface solver_map_eval
    procedure __APPEND(solver_map_eval_scalar,__PREC)
end interface

interface solver_map_eval
    procedure __APPEND(solver_map_eval_diag,__PREC)
end interface

interface solver_map_eval
    procedure __APPEND(solver_map_eval_matrix,__PREC)
end interface


interface solver_map_eval_inverse
    procedure __APPEND(solver_map_eval_inverse_scalar,__PREC)
end interface

interface solver_map_eval_inverse
    procedure __APPEND(solver_map_eval_inverse_diag,__PREC)
end interface

interface solver_map_eval_inverse
    procedure __APPEND(solver_map_eval_inverse_matrix,__PREC)
end interface

