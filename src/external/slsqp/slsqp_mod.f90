module slsqp_mod

    use, intrinsic :: iso_fortran_env
    
    use slsqp_mod_real64, only: slsqp_real64 => slsqp, slsqp_data_real64 => slsqp_data
    
    private
    public :: slsqp, slsqp_data_real64
    
    interface slsqp
        procedure :: slsqp_real64
    end interface

end module
