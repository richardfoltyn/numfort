      module slsqp_mod_real64
        use, intrinsic :: iso_fortran_env
        
        implicit none
        private
        
        public :: slsqp
        
        integer, parameter :: PREC = real64
        
        contains
        
        include "slsqp_optmz.f"
      end module
