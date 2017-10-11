module slsqpi_mod

    use, intrinsic :: iso_fortran_env
    
    use slsqp_mod_real64, only: slsqp_real64 => slsqp
    
    interface slsqp
        procedure :: slsqp_real64
    end interface

end module
