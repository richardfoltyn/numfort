

interface bernstein_check_input
    !*  Routines that perform input checks that are shared across all
    !   routines using polynonomials expressed in terms of the Bernstein basis.
    procedure __APPEND(bernstein_check_input,__PREC)
end interface

interface ppolyder
    procedure __APPEND(bernstein_ppolyder,__PREC)
end interface

interface ppolyder_impl
    procedure __APPEND(bernstein_ppolyder_impl,__PREC)
end interface

interface ppolyfit
   procedure __APPEND(bernstein_fit_deriv,__PREC)
end interface

interface bernstein_fit_deriv_check_input
   procedure __APPEND(bernstein_fit_deriv_check_input,__PREC)
end interface

interface bernstein_fit_deriv_impl
    procedure __APPEND(bernstein_fit_deriv_impl,__PREC)
end interface

interface ppolyval
    procedure __APPEND(bernstein_ppolyval_1d,__PREC)
end interface

interface ppolyval
    procedure __APPEND(bernstein_ppolyval_scalar,__PREC)
end interface

interface ppolyval_eval_impl
    procedure __APPEND(bernstein_ppolyval_eval_impl_scalar,__PREC)
end interface

interface ppolyval_eval_impl
    procedure __APPEND(bernstein_ppolyval_eval_impl_1d,__PREC)
end interface

interface ppolyval_eval_impl_deg2
    procedure __APPEND(bernstein_ppolyval_eval_impl_deg2_1d,__PREC)
end interface

interface ppolyval_eval_impl_deg3
    procedure __APPEND(bernstein_ppolyval_eval_impl_deg3_1d,__PREC)
end interface

interface ppolyval_eval_impl_degk
    procedure __APPEND(bernstein_ppolyval_eval_impl_degk_1d,__PREC)
end interface

interface bernstein_polyval
    procedure __APPEND(bernstein_polyval,__PREC)
end interface

interface ppoly_transform_basis
    procedure __APPEND(bernstein2power,__PREC)
end interface

!interface bernstein2power_deg3
!    procedure __APPEND(bernstein2power_deg3,__PREC)
!end interface
!
!interface bernstein2power_deg2
!    procedure __APPEND(bernstein2power_deg2,__PREC)
!end interface

