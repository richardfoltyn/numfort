real (PREC), intent(in), dimension(:) :: x
    !!  Array x where function (PDF, CDF, etc.) should be evaluated.
real (PREC), intent(out), dimension(:) :: fx
    !!  Array where function values are stored.
real (PREC), intent(in), optional :: loc
    !!  Location parameter (if applicable)
real (PREC), intent(in), optional :: scale
    !!  Scale parameter (if applicable)
real (PREC), intent(in), optional :: shape
    !! Shape parameter (if applicable)
