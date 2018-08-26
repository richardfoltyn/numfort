      abstract interface
        subroutine fobj_if (x, fx)
            import PREC
            real (PREC), intent(in), dimension(:), contiguous :: x
            real (PREC), intent(out), dimension(:), contiguous :: fx
        end subroutine
      end interface
