

module numfort_core_libc
    !*  Module implements wrappers around various math functions defined in
    !   the C standard library.

    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env

    implicit none
    private


    public :: log1p
    public :: expm1

    interface
        function c_log1p_double (value) result(res) bind(C, name='log1p')
            import
            real (C_DOUBLE), intent(in), value :: value
            real (C_DOUBLE) :: res
        end function

        function c_log1p_float (value) result(res) bind(C, name='log1pf')
           import
            real (C_FLOAT), intent(in), value :: value
            real (C_FLOAT) :: res
        end function
    end interface

    interface
        function c_expm1_double (value) result(res) bind(C, name='expm1')
            import
            real (C_DOUBLE), intent(in), value :: value
            real (C_DOUBLE) :: res
        end function

        function c_expm1_float (value) result(res) bind(C, name='expm1f')
            import
            real (C_FLOAT), intent(in), value :: value
            real (C_FLOAT) :: res
        end function

    end interface

    interface c_log1p
        procedure c_log1p_float, c_log1p_double
    end interface

    interface log1p
        procedure log1p_real32, log1p_real64
    end interface

    interface c_expm1
        procedure c_expm1_float, c_expm1_double
    end interface

    interface expm1
        procedure expm1_real32, expm1_real64
    end interface

    contains


impure elemental function log1p_real32 (x) result(res)
    !*  LOG1P implements the LOG(1+x) function
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: x
    real (PREC) :: res

    res = c_log1p (x)
end function

impure elemental function log1p_real64 (x) result(res)
    !*  LOG1P implements the LOG(1+x) function
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: x
    real (PREC) :: res

    res = c_log1p (x)
end function



impure elemental function expm1_real32 (x) result(res)
    !*  EXPM1 implements the EXP(x)-1 function
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: x
    real (PREC) :: res

    res = c_expm1 (x)
end function

impure elemental function expm1_real64 (x) result(res)
    !*  EXPM1 implements the EXP(x)-1 function
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: x
    real (PREC) :: res

    res = c_expm1 (x)
end function

end module
