MODULE simplex_csiro_real32

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(6, 30)
PRIVATE
PUBLIC :: minim, functn_if

INTERFACE
  SUBROUTINE functn_if(p, func)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(6, 30)
    REAL (dp), INTENT(IN)  :: p(:)
    REAL (dp), INTENT(OUT) :: func
  END SUBROUTINE
END INTERFACE

CONTAINS

    include "minim.f90"

END MODULE