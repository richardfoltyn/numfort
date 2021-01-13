

subroutine gaussian_kde (x, fhat, bw, trans_x)
    !*  Multivariate kernel densitity estimator with Gaussian kernel.
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  Array of points. By default, each row is assumed to contain one
        !   observation.
    real (PREC), intent(out), dimension(:), contiguous :: fhat
        !*  Kernel density estimate for each point in X.
    real (PREC), intent(in), optional :: bw
        !*  Optional bandwidth parameter. If not present, Scott's bandwidth
        !   rule is used.
    logical, intent(in), optional :: trans_x
        !*  If present and true, X is provided in transposed form, ie.
        !   each column contains one observation.

    logical :: ltrans_x
    real (PREC), dimension(:,:), allocatable :: x_T

    ltrans_x = .false.
    if (present(trans_x)) ltrans_x = trans_x

    if (trans_x) then
        allocate (x_T(size(x,2),size(x,1)))

        x_T(:,:) = transpose(x)

        call gaussian_kde_impl (x_T, fhat, bw)

        deallocate (x_T)
    else
        call gaussian_kde_impl (x, fhat, bw)
    end if

end subroutine



subroutine gaussian_kde_impl (x, fhat, bw)
    !*  Multivariate kernel densitity estimator with Gaussian kernel.
    !
    !   Implementation routine
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  Array of points. Each row is assumed to contain one observation.
    real (PREC), intent(out), dimension(:), contiguous :: fhat
    real (PREC), intent(in), optional :: bw

    real (PREC) :: fhat_i, x_ij, dist_scale, fhat_scale
    real (PREC), dimension(:), allocatable :: dist_i

    ! n = number of obs., d = number of vars
    integer :: n, d, i, j, k
    ! local bandwidth
    real (PREC) :: lbw

    n = size(x, 1)
    d = size(x, 2)

    ! Scott's bandwidth rule, used as default
    lbw = n ** (-1.0_PREC/(d + 4.0_PREC))
    if (present(bw)) lbw = bw

    ! Scaling factor for squared distance
    dist_scale = 1.0_PREC / (2.0_PREC * lbw ** 2.0_PREC)
    fhat_scale =  1.0_PREC / (n * (2 * PI) ** (d/2.0_PREC) * lbw ** d)

    !$omp parallel default(none) &
    !$omp shared(x,fhat,fhat_scale,dist_scale,n,d) &
    !$omp private(i,j,k,fhat_i,dist_i,x_ij)

    allocate (dist_i(n))

    !$omp do schedule(static)
    do i = 1, n
        dist_i(:) = 0.0
        ! Compute squared distance to any other point indexed by k;
        ! For efficient memory access sum over columns first.
        do j = 1, d
            x_ij = x(i,j)
            do k = 1, n
                dist_i(k) = dist_i(k) + (x(k,j) - x_ij) ** 2.0_PREC
            end do
        end do

        ! Compute density estimate at point x_i
        fhat_i = 0.0
        do k = 1, n
            fhat_i = fhat_i + exp( - dist_i(k) * dist_scale)
        end do

        fhat(i) = fhat_i * fhat_scale
    end do
    !$omp end do

    deallocate (dist_i)

    !$omp end parallel

end subroutine

