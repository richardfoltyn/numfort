      abstract interface
        subroutine fobj_if (n, m, x, fx)
            import PREC
            integer, intent(in) :: n
            integer, intent(in) :: m
            real (PREC), intent(in), dimension(*) :: x
            real (PREC), intent(out), dimension(*) :: fx
        end subroutine
      end interface
