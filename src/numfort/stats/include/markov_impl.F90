


pure function __APPEND(is_trans_matrix,__PREC) (mat, tol, transposed) result(res)
    !*  IS_TRANS_MATRIX checks whether argument is a valid transition matrix
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: mat
        !*  Transition matrix to check
    real (PREC), intent(in), optional :: tol
        !*  Tolerance on row sums (default: 1e-12)
    logical, intent(in), optional :: transposed
        !*  If true, assume that transition matrix is transposed, ie. columns
        !   sum to 1.0 instead of rows.
    logical :: res

    real (PREC) :: ltol
    integer :: dim
    logical :: ltransposed

    ltol = 1e-12_PREC
    ltransposed = .false.

    if (present(tol)) ltol = tol
    if (present(transposed)) ltransposed = transposed

    ! Dimension across which elements should sum to 1.0
    dim = 2
    if (ltransposed) dim = 1

    ! Check that all elements are non-negative and that rows (or columns if
    ! transposed=.true.) approximately sum to 1.0
    res = all(mat >= 0.0_PREC) .and. all(abs(sum(mat, dim=dim) - 1.0_PREC) < ltol)

end function


subroutine __APPEND(markov_approx_input_checks,__PREC) (rho, sigma, states, &
        tm, status)
    !*  MARKOV_APPROX_INPUT_CHECKS validates inputs for the ROUWENHORST
    !   and TAUCHEN routines.

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in) :: rho
    real (PREC), intent(in) :: sigma
    real (PREC), intent(in), dimension(:) :: states
    real (PREC), intent(in), dimension(:,:) :: tm
    type (status_t), intent(out) :: status

    integer :: n

    status = NF_STATUS_OK
    n = size(states)

    ! Input checks
    if ((abs(rho) >= 1.0_PREC) .and. (n > 1)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if ((sigma < 0.0_PREC) .and. (n > 1)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (n /= size(tm,1) .or. n /= size(tm,2)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

100 continue

end subroutine

subroutine __APPEND(rouwenhorst,__PREC) (rho, sigma, states, tm, sigma_cond, status)
    !*  ROUWENHORST returns the discretized approximation of an AR(1) process
    !   using the Rouwenhorst method.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in) :: rho
        !*  Autocorrelation coefficient
    real (PREC), intent(in) :: sigma
        !*  Conditional or unconditional standard deviation
        !   (depending on SIGMA_COND)
    real (PREC), intent(out), dimension(:) :: states
        !*  Contains discretized state vector on successful exit.
    real (PREC), intent(out), dimension(:,:) :: tm
        !*  Contains transition matrix on successful exit.
    logical, intent(in), optional :: sigma_cond
        !*  If present and true, argument SIGMA is interpreted as the
        !   conditional standard deviation (default: true).
    type (status_t), intent(out), optional :: status
        !*  Exit status (optional)

    type (status_t) :: lstatus
    integer :: n, i, j, m
    real (PREC), dimension(:,:), allocatable :: work
    real (PREC) :: eps, p, sigma_z
    logical :: lsigma_cond

    lstatus = NF_STATUS_OK

    n = size(states)

    ! Check inputs
    call markov_approx_input_checks (rho, sigma, states, tm, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Permit calling routine with size-0 arrays, do nothing in that case
    if (n == 0) goto 100

    lsigma_cond = .true.
    if (present(sigma_cond)) lsigma_cond = sigma_cond

    ! Treat degenerate case with one state separately
    if (n == 1) then
        tm(1,1) = 1.0_PREC
        states(1) = 0.0_PREC
        goto 100
    end if

    if (lsigma_cond) then
        sigma_z = sigma / sqrt(1.0_PREC - rho ** 2.0_PREC)
    else
        sigma_z = sigma
    end if

    ! Construct discretized state space
    eps = sqrt(n-1.0_PREC) * sigma_z
    call linspace (states, -eps, eps)

    ! At this point the state space has at least N=2
    allocate (work(n+1,n+1), source=0.0_PREC)

    p = (1.0_PREC+rho)/2.0_PREC

    ! Initial 2-by-2 transition matrix
    tm(1,1) = p
    tm(1,2) = 1.0_PREC - p
    tm(2,1) = 1.0_PREC - p
    tm(2,2) = p

    do m = 2, n-1
        call rouwenhorst_pad_matrix (tm(1:m,1:m), work)

        do j = 1, m+1
            do i = 1, m+1
                tm(i,j) = p * work(1+i,1+j) &
                    + (1.0_PREC - p) * work(1+i,j) &
                    + (1.0_PREC - p) * work(i,j+1) &
                    + p * work(i,j)
            end do
        end do

        do j = 1, m+1
            tm(2:m,j) = tm(2:m,j) / 2.0_PREC
        end do
    end do

    ! normalize row sums to 1
    do i = 1, n
        tm(i,:) = tm(i,:) / sum(tm(i,:))
    end do

100 continue
    if (present(status)) status = lstatus

    if (allocated(work)) deallocate (work)

end subroutine


subroutine __APPEND(rouwenhorst_pad_matrix,__PREC) (mat, res)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: mat
    real (PREC), intent(out), dimension(:,:) :: res

    integer :: n, i

    n = size(mat, 1)

    res(1:n+2,1:n+2) = 0.0_PREC
    do i = 1, n
        res(2:n+1,i+1) = mat(1:n,i)
    end do
end subroutine



subroutine __APPEND(tauchen,__PREC) (rho, sigma, states, tm, m, sigma_cond, &
        status)
    !*  TAUCHEN returns the discretized approximation of an AR(1) process
    !   using the Tauche (1986) method.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: rho
        !*  Autocorrelation coefficient
    real (PREC), intent(in) :: sigma
        !*  Conditional (default) or unconditional standard deviation
        !   (interpretation depends on the argument SIGMA_COND)
    real (PREC), intent(out), dimension(:) :: states
        !*  Contains discretized state space on successful exit.
    real (PREC), intent(out), dimension(:,:) :: tm
        !*  Contains transition matrix on successful exit.
    real (PREC), intent(in), optional :: m
        !*  If present, determines the range of the discretized state space
        !   as (sigma * M), where SIGMA is the unconditional standard deviation
        !   (default: 3).
    logical, intent(in), optional :: sigma_cond
        !*  If present and true, interpret SIGMA argument as the conditional
        !   standard deviation (default: true).
    type (status_t), intent(out), optional :: status
        !*  Exit status (optional)

    logical :: lsigma_cond
    real (PREC) :: lm
    real (PREC) :: sigma_z, sigma_e
    real (PREC) :: zjk, w2, zj1, zjn
    integer :: i, k, n
    type (status_t) :: lstatus

    type (__APPEND(dnorm,__PREC)), parameter :: norm = &
        __APPEND(dnorm,__PREC) (loc=0.0_PREC,scale=1.0_PREC)

    lstatus = NF_STATUS_OK

    ! Input checks
    call markov_approx_input_checks (rho, sigma, states, tm, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    n =  size(states)

    if (n == 0) then
        ! Permit calling routine with size-0 arrays; do nothing in that case
        goto 100
    else if (n == 1) then
        ! Degenerate process with single state and certain transition
        states(1) = 0.0
        tm(1,1) = 1.0
        goto 100
    end if

    lsigma_cond = .true.
    lm = 3

    if (present(sigma_cond)) lsigma_cond = sigma_cond
    if (present(m)) lm = m

    if (lsigma_cond) then
        sigma_z = sigma / sqrt(1.0_PREC - rho ** 2.0_PREC)
        sigma_e = sigma
    else
        sigma_z = sigma
        sigma_e = sigma_z * sqrt(1.0_PREC - rho ** 2.0_PREC)
    end if

    call linspace(states, -sigma_z * lm, sigma_z * lm)
    w2 = (states(2) - states(1)) / 2.0_PREC

    do i = 1, n
        do k = 2, n-1
            zjk = states(k) - rho * states(i)

            tm(i, k) = cdf (norm, (zjk + w2) / sigma_e) &
                - cdf (norm, (zjk - w2) / sigma_e)
        end do

        zj1 = states(1) - rho * states(i)
        zjn = states(n) - rho * states(i)

        tm(i, 1) = cdf (norm, (zj1 + w2) / sigma_e)
        tm(i, n) = 1 - cdf (norm, (zjn - w2) / sigma_e)
    end do

    ! normalize row sums to 1
    do i = 1, n
        tm(i,:) = tm(i,:) / sum(tm(i,:))
    end do

100 continue
    if (present(status)) status = lstatus

end subroutine


subroutine __APPEND(ergodic_dist,__PREC)  (tm, edist, inverse, maxiter, &
        is_transposed, status)
    !*  ERGODIC_DIST returns the ergodic distribution implied by a given
    !   Markov transition matrix using one of two methods (see argument INVERSE).
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:,:) :: tm
        !*  Markov transition matrix.
    real (PREC), intent(out), dimension(:) :: edist
        !*  Contains ergodic distribution on successful exit
    logical, intent(in), optional :: inverse
        !*  If present and true, compute ergodic distribution using the
        !   "inversion" method (default: true). Otherwise, iterate on the
        !   transition matrix until element-wise convergence.
    integer, intent(in), optional :: maxiter
        !*  Maximum number of iterations if INVERSE is .FALSE.
    logical, intent(in), optional :: is_transposed
        !*  If present and true, assume that TM is a transposed transition
        !   matrix, ie TM(i,j) = Prob(x'=i|x=j). (default: false)
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    ! Working arrays
    real (PREC), dimension(:,:), allocatable :: tm_inv, tm_T
    real (PREC), dimension(:), allocatable, target :: mu1, mu2
    real (PREC), dimension(:), pointer, contiguous :: ptr1, ptr2, ptr3

    integer :: i, lmaxiter, n
    logical :: linverse, lis_transposed
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    n = size(tm, 1)

    if (size(edist) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    linverse = .true.
    lmaxiter = 10000
    lis_transposed = .false.
    if (present(inverse)) linverse = inverse
    if (present(maxiter)) lmaxiter = maxiter
    if (present(is_transposed)) lis_transposed = is_transposed

    if (lis_transposed) then
        allocate (tm_T(n,n), source=tm)
    else
        allocate (tm_T(n, n))
        tm_T(:,:) = transpose(tm)
    end if

    ! compute ergodic distribution using "inverse" method
    if (linverse) then

        forall (i=1:n) tm_T(i,i) = tm_T(i,i) - 1

        tm_T(n, :) = 1

        allocate (tm_inv(n,n))

        call inv(tm_T, tm_inv, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        edist = tm_inv(:, n)
    else
        ! compute ergodic distribution by iteration
        allocate (mu1(n), mu2(n))
        ! initialize with uniform distribution over states
        mu1 = 1.0_PREC / n

        ptr1 => mu1
        ptr2 => mu2

        do i = 1, lmaxiter
            ptr2 = matmul(tm_T, ptr1)

            ! check whether convergence in mu was achieved
            if (all(abs(ptr1 - ptr2) < 1d-12)) exit

            ptr3 => ptr1
            ptr1 => ptr2
            ptr2 => ptr3
        end do

        edist = ptr2
    end if

    ! normalize to 1
    edist = edist / sum(edist)

    if (abs(sum(edist) - 1.0_PREC) > 1e-12_PREC) then
        lstatus = NF_STATUS_INVALID_STATE
        goto 100
    end if

    if (any(edist < 0.0_PREC)) then
        lstatus = NF_STATUS_INVALID_STATE
        goto 100
    end if

100 continue

    if (present(status)) status = lstatus

    if (allocated(tm_T)) deallocate (tm_T)
    if (allocated(tm_inv)) deallocate (tm_inv)
    if (allocated(mu1)) deallocate (mu1)
    if (allocated(mu2)) deallocate (mu2)

end subroutine


subroutine __APPEND(moments,__PREC) (states, tm, mean, acorr, sigma, sigma_e, &
        edist, inverse, status)
    !*  MOMENTS computes selected moments implied by a Markov
    !   process with given state space and transition matrix.

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: states
        !*  Markov chain state space
    real (PREC), intent(in), dimension(:,:) :: tm
        !*  Transition matrix with element (i,j) being the transition
        !   probability Prob[x_t=s_j|x_{t-1}=s_i]
    real (PREC), intent(out), optional :: mean
        !*  If present, contains unconditional mean on successful exit.
    real (PREC), intent(out), optional :: acorr
        !*  On exit, contains implied autocorrelation coefficient
    real (PREC), intent(out), optional :: sigma
        !*  On exit, contains implied unconditional std. dev.
    real (PREC), intent(out), optional :: sigma_e
        !*  On exit, contains implied std. dev. of AR(1) error
        !   term (ie. conditional std. dev.)
    real (PREC), intent(in), dimension(:), optional, target :: edist
        !*  If present, contains the ergodic distribution implied by transition
        !   matrix. Will be computed if missing.
    logical, intent(in), optional :: inverse
        !*  Use inverse of transition matrix when computing ergodic distribution
        !   (ignored if EDIST argument is present).
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:), pointer :: ptr_pmf
    real (PREC), dimension(:), allocatable :: x
    logical :: linverse
    real (PREC) :: mean_x, var_x, covar_x, acorr_x
    integer :: i, j, n
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    n = size(states)

    linverse = .true.
    if (present(inverse)) linverse = inverse

    call assert_alloc_ptr (edist, n, ptr_pmf)

    if (.not. present(edist)) then
        call markov_ergodic_dist (tm, ptr_pmf, inverse=linverse, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! Unconditional mean
    mean_x = dot_product (ptr_pmf, states)

    ! Demeaned states
    allocate (x(n), source=states)
    x(:) = x - mean_x

    ! Unconditional variance
    var_x = 0.0
    do i = 1, n
        var_x = var_x + ptr_pmf(i) * x(i) ** 2
    end do

    ! Compute autocovariane
    covar_x = 0.0
    do j = 1, n
        do i = 1, n
            covar_x = covar_x + x(i) * x(j) * tm(i,j) * ptr_pmf(i)
        end do
    end do

    if (var_x > 0.0_PREC) then
        acorr_x = covar_x / var_x
    else
        ! If var_x == 0 we have a degenerate distribution with only one
        ! state, hence there is perfect correlation.
        acorr_x = 1.0_PREC
    end if

    if (present(mean)) mean = mean_x
    if (present(sigma_e)) sigma_e = sqrt((1.0_PREC-acorr_x**2.0_PREC) * var_x)
    if (present(sigma)) sigma = sqrt(var_x)
    if (present(acorr)) acorr = acorr_x

100 continue

    if (present(status)) status = lstatus

    call assert_dealloc_ptr (edist, ptr_pmf)
    if (allocated(x)) deallocate(x)

end subroutine
