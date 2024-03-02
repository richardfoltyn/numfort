

subroutine gini_impl (x, pmf, g, assume_sorted)
    !*  GINI returns the Gini coefficient for a discrete probabability distribution
    !   over states X with given PMF.
    !
    !   The formula used follows Wikipedia
    !   https://en.wikipedia.org/wiki/Gini_coefficient#Discrete_probability_distribution
    !
    real (PREC), intent(in), dimension(:), contiguous, target :: x
    real (PREC), intent(in), dimension(:), contiguous, target :: pmf
    real (PREC), intent(out) :: g
    logical, intent(in), optional :: assume_sorted

    integer, dimension(:), allocatable :: iorder
    real (PREC), dimension(:), allocatable, target :: x_sort, pmf_sort
    real (PREC), dimension(:), contiguous, pointer :: ptr_x, ptr_pmf
    logical :: lassume_sorted

    real (PREC) :: s1, s2
    integer :: n, i
    lassume_sorted = .false.


    if (present(assume_sorted)) lassume_sorted = assume_sorted

    g = ieee_value (0.0_PREC, IEEE_SIGNALING_NAN)

    n = size (x)
    if (size(x) /= size (pmf)) return

    if (n == 1) then
        ! For degenerate distributions the Gini must be 0.
        ! Take a short-cut as this also avoids issues if size(x) = 1 and x(1) = 0,
        ! as then the formula below returns NaN.
        g = 0.0
        return
    end if

    ! --- Ensure sorted inputs ---

    if (assume_sorted) then
        ! Arrays are already sorted, point to dummy arguments
        ptr_x => x
        ptr_pmf => pmf
    else
        ! Sort arrays
        allocate (iorder(n), x_sort(n), pmf_sort(n))

        ! Get sort order
        call argsort (x, iorder)

        do i = 1, n
            x_sort(i) = x(iorder(i))
            pmf_sort(i) = pmf(iorder(i))
        end do
        
        ptr_x => x_sort
        ptr_pmf => pmf_sort
    end if

    ! --- Compute Gini ---

    g = 0.0
    s1 = 0.0
    s2 = 0.0

    do i = 1, n
        s2 = s1 + ptr_x(i) * ptr_pmf(i)
        g = g + (s1 + s2) * ptr_pmf(i)
        s1 = s2
    end do

    g = 1.0_PREC - g / s2

end subroutine