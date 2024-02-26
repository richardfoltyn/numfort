



pure subroutine h12_compute (lpivot, l, u, h)
    integer, intent(in) :: lpivot
    integer, intent(in) :: l
    real (PREC), intent(inout), dimension(:) :: u
    real (PREC), intent(out) :: h

    integer :: m, i
    real (PREC) :: vp, s, umax

    m = size(u)

    ! Input checks: we need to enforce that
    !   1 <= lpivot < l1 <= m
    ! Don't report any error since in the original version this was
    ! used to effectively apply an identity transformation
    if ((lpivot < 1) .or. (lpivot >= l) .or. (l > m)) return

    ! Step 11(b)
    ! Underflow-resistant computation of sum of squares
    umax = abs(u(lpivot))
    do i = l, m
        umax = max(umax, abs(u(i)))
    end do

    if (umax <= 0.0_PREC) return

    umax = 1.0_PREC / umax

    vp = u(lpivot)

    ! Step 1:
    s = (vp * umax)**2.0_PREC
    do i = l, m
        s = s + (u(i) * umax)**2.0_PREC
    end do
    s = sqrt(s) / umax

    ! Step 2: if vp > 0, s = -s
    if (vp > 0.0_PREC) s = -s

    ! Step 3: h = vp - s, vp = s
    h = vp - s
    u(lpivot) = s

end subroutine



pure subroutine h12_apply_1d (lpivot, l, u, h, c, trans, status, msg)
    integer, intent(in) :: lpivot
    integer, intent(in) :: l
    real (PREC), intent(inout), dimension(:) :: u
    real (PREC), intent(in) :: h
    real (PREC), intent(out), dimension(:), target :: c
    logical, intent(in), optional :: trans
        !*  Ignored. Only present for API-compatibility with implementation
        !   for 2d arrays.
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    logical :: ltrans
    real (PREC), dimension(:,:), pointer :: ptr_c

    ltrans = .false.
    if (present(trans)) ltrans = trans

    if (ltrans) then
        ptr_c(1:1,1:size(c)) => c
    else
        ptr_c(1:size(c), 1:1) => c
    end if

    call h12_apply (lpivot, l, u, h, ptr_c, ltrans, status, msg)

end subroutine



pure subroutine h12_apply_2d (lpivot, l, u, h, c, trans, status, msg)
    integer, intent(in) :: lpivot
    integer, intent(in) :: l
    real (PREC), intent(inout), dimension(:) :: u
    real (PREC), intent(in) :: h
    real (PREC), intent(out), dimension(:,:) :: c
    logical, intent(in), optional :: trans
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    real (PREC) :: b, vp, s
    integer :: m, i, j
    logical :: ltrans
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    ltrans = .false.
    if (present(trans)) ltrans = trans

    m = size(u)

    ! Input checks: we need to enforce that
    !   1 <= lpivot < l1 <= m
    ! Don't report any error since in the original version this was
    ! used to effectively apply an identity transformation
    if ((lpivot < 1) .or. (lpivot >= l) .or. (l > m)) goto 100

    ! Quietly exit if there is nothing to do
    if (size(c) == 0) goto 100

    if ((ltrans .and. m /= size(c, 2)) .or. &
            (.not. ltrans .and. m /= size(c, 1))) then
        lstatus = NF_STATUS_INVALID_ARG
        if (present(msg)) then
            msg = "H12_APPLY: Non-conformable input arrays A, C"
        end if
        goto 100
    end if

    ! Apply the transformation I + U*(U^T)/B to C
    ! Step 5: Set b := vp * h
    vp = u(lpivot)
    b = vp * h

    ! B must be nonpositive here. If B = 0, return
    if (b >= 0.0_PREC) goto 100

    ! Step 7: Apply transformation to C

    ! Step 8
    if (ltrans) then
        do j = 1, size(c, 1)
            s = c(j,lpivot) * h
            do i = l, m
                s = s + c(j,i)*u(i)
            end do
            s = s / b

            ! Step 9: set c(p,j) = c(p,j) + s * h
            c(j,lpivot) = c(j,lpivot) + s * h

            ! Step 10: Set c(i,j) = c(i,j) + s * u(i) for i = l,...,m
            do i = l, m
                c(j,i) = c(j,i) + s * u(i)
            end do
        end do
    else
        do j = 1, size(c, 2)
            s = c(lpivot,j) * h
            do i = l, m
                s = s + c(i,j)*u(i)
            end do
            s = s / b

            ! Step 9: set c(p,j) = c(p,j) + s * h
            c(lpivot,j) = c(lpivot,j) + s * h

            ! Step 10: Set c(i,j) = c(i,j) + s * u(i) for i = l,...,m
            do i = l, m
                c(i,j) = c(i,j) + s * u(i)
            end do
        end do
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine h12_1d (imode, lpivot, l, u, h, c, trans, status, msg)
    integer, intent(in) :: imode
    integer, intent(in) :: lpivot
    integer, intent(in) :: l
    real (PREC), intent(inout), dimension(:) :: u
    real (PREC), intent(out) :: h
    real (PREC), intent(out), dimension(:) :: c
    logical, intent(in), optional :: trans
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    if (present(status)) status = NF_STATUS_OK

    if (imode == 1) then
        call h12_compute (lpivot, l, u, h)
    end if

    call h12_apply (lpivot, l, u, h, c, trans, status, msg)

end subroutine



pure subroutine h12_2d (imode, lpivot, l, u, h, c, trans, status, msg)
    integer, intent(in) :: imode
    integer, intent(in) :: lpivot
    integer, intent(in) :: l
    real (PREC), intent(inout), dimension(:) :: u
    real (PREC), intent(out) :: h
    real (PREC), intent(out), dimension(:,:), optional :: c
    logical, intent(in), optional :: trans
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    if (present(status)) status = NF_STATUS_OK

    if (imode == 1) then
        call h12_compute (lpivot, l, u, h)
    end if

    if (present(c)) then
        call h12_apply (lpivot, l, u, h, c, trans, status, msg)
    end if

end subroutine



subroutine hfti_1d (a, b, tau, krank, rnorm, w, ip, status, msg)
    !*  HFTI is the wrapper around HFTI for 1-dimensional arguments B.
    real (PREC), intent(inout), dimension(:,:), contiguous :: a
    real (PREC), intent(inout), dimension(:), contiguous, target :: b
    real (PREC), intent(in) :: tau
    integer, intent(out) :: krank
    real (PREC), intent(out) :: rnorm
    real (PREC), intent(out), dimension(:), contiguous :: w
    integer, intent(inout), dimension(:), contiguous :: ip
    type (status_t), intent(inout), optional :: status
    character (*), intent(out), optional :: msg

    real (PREC) :: rnorm1d(1)
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_b

    ! Ensure that value does not change if routine exists prematurely
    rnorm1d(1) = rnorm

    ptr_b(1:size(b),1:1) => b

    call hfti (a, ptr_b, tau, krank, rnorm1d, w, ip, status, msg)

    rnorm = rnorm1d(1)

end subroutine



subroutine hfti_query (m, n, w, ip)
    integer, intent(in) :: m, n
    integer, intent(out) :: w, ip

    w = min(m, n) + n
    ip =  min(m, n)

end subroutine



subroutine hfti_2d (a, b, tau, krank, rnorm, w, ip, status, msg)
    !*  HFTI solves a linear least squares problem Ax = b
    !   by Householder transformations.
    !
    !   Implementation of algorithm 14.9 in
    !       Lawson, Hanson (1995): Solving Least Squares Problems
    !   Based on the Fortran code provided on NETLIB for the 1995 reprint.
    !
    !   Implementation notes:
    !       -   No array dimensions are passed. These are inferred directly
    !           from the array arguments, input arrays should thus be
    !           subscripted as needed.
    !       -   Routine and in particular the routine H12 called from HTFI
    !           do not perform contiguous memory access, hence the
    !           CONTIGUOUS attribute must not be specified.
    real (PREC), intent(inout), dimension(:,:), contiguous :: a
        !*  Coefficient matrix A initially contains the M x N matrix
        !   of the least squares problem Ax = b. Both M >= N and
        !   N < M are permitted. Contents of A will be modified by
        !   the subroutine.
    real (PREC), intent(inout), dimension(:,:), contiguous :: b
        !*  Array b of least-squares problem Ax = b. Must be of size
        !   of at least max(M,N). On exit, b will contain
        !   the vector x if length N (or vectors, if size(b, 2) > 1).
    real (PREC), intent(in) :: tau
        !*  Absolute tolerance parameter for pseudorank determination.
    integer, intent(out) :: krank
        !*  Contains pseudorank of A on exit.
    real (PREC), intent(out), dimension(:), contiguous :: rnorm
        !*  On exit, contains the Euclidean norm of the residual vector
        !   for the problem defined by the corresponding column of b.
    real (PREC), intent(out), dimension(:), contiguous, target :: w
        !*  Working array. Combines working arrays H and G from original
        !   implementation into one, such that W = [H, G].
    integer, intent(inout), dimension(:), contiguous :: ip
        !*  Array containing indices describing the permutation of column
        !   vectors.
    type (status_t), intent(out), optional :: status
        !*  Optional exit status flag.
    character (*), intent(out), optional :: msg
        !*  Optional status message. Changed by routine only on
        !   unsuccessful termination.

    integer :: m, n, mb, nb, k, mu, lmax, i, j, l, jb
    real (PREC) :: hmax, diff, tmp, s
    real (PREC), parameter :: FACTOR = 0.001_PREC
    type (status_t) :: lstatus
    real (PREC), dimension(:), pointer, contiguous :: g, h
        ! Pointers to arrays G, H as used in original code

    lstatus = NF_STATUS_OK

    ! Initialize to avoid uninitialized variable compiler warnings
    krank = -1

    m = size(a, 1)
    n = size(a, 2)
    mb = size(b, 1)
    nb = size(b, 2)
    mu = min(m, n)

    ! === Check inputs ===

    if (mu <= 0) then
        ! Special case: don't flag 0-size input arrays as error to
        ! be compatible with original implementation.
        krank = 0
        goto 100
    end if

    ! B will be used to store x from Ax = b, so it needs to have
    ! dimension max(M,N) since x has dimension N
    call assert_input (max(m,n) == mb, 'size(B,1) = max(m,n) required')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(rnorm) == nb, 'size(RNORM) = size(B,2) required')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(ip) >= mu, 'size(IP) >= min(M,N) required')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(w) >= (mu + n), 'size(W) >= (min(M,N) + N) required')
    if (lstatus /= NF_STATUS_OK) goto 100

    ! === Algorith HFTI ===

    krank = 0
    hmax = 0.0_PREC
    rnorm = 0.0
    diff = 0.0

    ! Associate points to replace arrays H, G in original code
    h(1:mu) => w(1:mu)
    g(1:n) => w(mu+1:mu+n)

    do j = 1, mu

        ! === Determine LMAX ===

        ! This combines steps 4-8 in the algorithm.

        ! StepS 4+5: Update squared column lengths and find LMAX
        if (j > 1) then
            lmax = j
            do l = j, n
                h(l) = h(l) - a(j-1,l)**2.0_PREC
                if (h(l) > h(lmax)) lmax = l
            end do

            ! Step 6: Check difference
            diff = (hmax + factor*h(lmax)) - hmax
        end if

        ! Step 7+8: Executed for J = 1 or if DIFF = 0 was found above
        ! Compute squared column lengths and find LMAX
        if (j == 1 .or. diff <= 0.0_PREC) then
            lmax = j
            do l = j, n
                h(l) = a(j,l) ** 2.0_PREC
                do i = j+1, m
                    h(l) = h(l) + a(i,l)**2.0_PREC
                end do
                if (h(l) > h(lmax)) lmax = l
            end do
            hmax = h(lmax)
        end if

        ! === LMAX determined ===

        ! Step 9
        ip(j) = lmax

        ! Step 10: Interchange columns J and LMAX of A, set H(LMAX) = H(J)
        if (ip(j) /= j) then
            do i = 1, m
                tmp = a(i,j)
                a(i,j) = a(i,lmax)
                a(i,lmax) = tmp
            end do

            ! Set H(LMAX) = H(J)
            h(lmax) = h(j)
        end if

        ! Step 11: Compute the j-th Householder transformation and apply to A
        call h12 (1, j, j+1, a(:,j), h(j), a(:,j+1:n), status=lstatus, msg=msg)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! Step 12: Apply Householder transformation to B
        call h12 (2, j, j+1, a(:,j), h(j), b(1:m,1:nb), status=lstatus, msg=msg)
        if (lstatus /= NF_STATUS_OK) goto 100

    end do

    ! Step 13: Determine pseudo-rank
    ! Diagonal elements of R, stored in A(1,1) through A(mu,mu) are
    ! non-increasing in magnitude, so we need not check beyond the
    ! first diagonal element A(j,j) <= TAU.
    k = mu
    do j = 1, mu
        if (abs(a(j,j)) <= tau) then
            k = j - 1
            exit
        end if
    end do

    ! Compute the norms of residual vectors
    do jb = 1, nb
        if (k+1 <= m) then
            rnorm(jb) = NRM2 (m-k, b(k+1:m,jb), 1)
        end if
    end do

    ! Handle special case of pseudo-rank = 0
    if (k == 0) then
        b = 0.0_PREC
        krank = 0
        goto 100
    end if


    ! Step 15: For K < N, determine orthogonal transformations K_i
    if (k < n) then
        ! Step 16: If K < N, compute Householder decomposition of the first
        ! K rows.
        do i = k, 1, -1
            call h12 (1, i, k+1, a(i,:), g(i), a(1:i-1,:), trans=.true., &
                status=lstatus, msg=msg)
            if (lstatus /= NF_STATUS_OK) goto 100
        end do
    end if

    do jb = 1, nb
        ! Step 17: x_k = B_k / A_kk
        b(k,jb) = b(k,jb) / a(k,k)

        ! Step 18
        do i = k-1, 1, -1
            s = 0.0_PREC
            do j = i + 1, k
                s = s + a(i,j)*b(j,jb)
            end do
            b(i,jb) = (b(i,jb) - s) / a(i,i)
        end do

        ! Step 19
        if (k < n) then

            ! Complete computation of solution vector
            ! Step 20
            b(k+1:n,jb) = 0.0_PREC

            do i = 1, k
                call h12 (2, i, k+1, a(i,:), g(i), b(1:n,jb), trans=.false., &
                    status=lstatus, msg=msg)
                if (lstatus /= NF_STATUS_OK) goto 100
            end do
        end if

        ! Step 22: Re-order the solution vector to compensate the column
        ! interchanges.
        do j = mu, 1, -1
            if (ip(j) == j) cycle

            l = ip(j)
            tmp = b(l,jb)
            b(l,jb) = b(j,jb)
            b(j,jb) = tmp
        end do
    end do

    krank = k

100 continue

    ! Write back ISTATUS if present
    if (present(status)) status = lstatus

    contains

    subroutine assert_input (cond, lmsg)
        logical, intent(in) :: cond
        character (*), intent(in) :: lmsg

        lstatus = NF_STATUS_OK
        if (.not. cond) then
            lstatus = NF_STATUS_INVALID_ARG

            if (present(msg)) then
                msg = 'HFTI: Invalid input: ' // lmsg
            end if
        end if
    end subroutine

end subroutine



subroutine nnls_query (m, n, w, z, index)
    integer, intent(in) :: m, n
    integer, intent(out) :: w, z, index

    w = n
    z = m
    index = n

end subroutine



subroutine nnls (a, b, x, rnorm, w, z, index, status, msg)
    !*  NNLS implements the non-negative least squares solver
    !   from Algorithm 23.10 in
    !       Lawson, Hanson (1995): Solving least squares problems
    !   which finds x such that
    !       Ax = b      s.t. ||x|| >= 0
    !
    real (PREC), intent(inout), dimension(:,:), contiguous :: a
        !*  Coefficient matrix A with dimensions (M,N). Either M >= N
        !   or M < N is permissible. There is no restriction on the rank
        !   of A.
        !   On termination, A contains the matrix A.Q, where Q
        !   is an orthogonal matrix implicitly generated by
        !   this routine.
    real (PREC), intent(inout), dimension(:), contiguous :: b
        !*  B initially contains the M-vector b. On termination, B
        !   contains Q.b.
    real (PREC), intent(out), dimension(:), contiguous :: x
        !*  On termination with IMODE = 1, X contains the solution
        !   vector.
    real (PREC), intent(out) :: rnorm
        !*  On termination, RNORM contains the Euclidean norm of the
        !   final residual vector RNORM = ||b - Ax||.
    real (PREC), intent(out), dimension(:), contiguous :: w
        !*  Array of at least size max(M,N).
        !   On termination, W contains the N-dimensional
        !   dual vector W = A^T (b-Ax)
    real (PREC), intent(out), dimension(:), contiguous :: z
        !*  Working space of at least size M.
    integer, intent(out), dimension(:), contiguous :: index
        !*  Integer working space of at least size N.
    type (status_t), intent(out), optional :: status
        !*  Optional exit status flag.
    character (*), intent(out), optional :: msg
        !*  Optional status message. Changed by routine only on
        !   unsuccessful termination.

    real (PREC), parameter :: FACTOR = 0.01_PREC
    integer :: m, n, np
    integer :: iter, i, iz1, iz2, iz, izmax, j, jj, jz, ip
    real (PREC) :: sm, wmax, tmp, up, diff, unorm, ztest, alpha, t
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    ! Original success exit code
    lstatus%code_orig = STATUS_OK

    ! Initialize to avoid uninitialized variable compiler warnings
    rnorm = huge(0.0_PREC)

    m = size(a, 1)
    n = size(a, 2)

    ! === Input checks ====

    call assert_input (size(a) > 0, 'size(A) > 0 required')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(b) == m, 'Non-conformable arrays A, B')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(w) >= n, 'size(W) >= N required')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(z) >= m, 'size(Z) >= M required')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(index) >= n, 'size(INDEX) >= N required')
    if (lstatus /= NF_STATUS_OK) goto 100

    ! === Solve least squares problem ===

    ! Step 1: initialize arrays INDEX and X
    x = 0.0_PREC
    do i = 1, n
        index(i) = i
    end do

    ! === Main loop ===

    ! Lower and upper bound of set Z
    iz1 = 1
    iz2 = N
    ! Size of set P
    np = 0
    ztest = 0.0

    ! Step 2: Compute coefficients of the dual (negative gradient) vector W
    do iz = iz1, iz2
        j = index(iz)
        w(j) = dot_product (a(np+1:m,j), b(np+1:m))
    end do

    ! Step 3: Termination condition if either of these conditions is met:
    !   (1) The set Z is empty, indicated by iz1 > iz2
    !   (2) w(j) <= 0 for all j in Z, indicated by NSETP >= M
    ! Quit if all coefficients are already in the solution, or if M columns
    ! of A have been triangularized
    main: do while (iz1 <= iz2 .and. np < m)

        ! Step 4: Find largest positive W(j)
        wmax = 0.0_PREC
        izmax = 1
        do iz = iz1, iz2
            j = index(iz)
            if (w(j) > wmax) then
                wmax = w(j)
                izmax = iz
            end if
        end do

        ! if WMAX <= 0.0, terminate loop. This indicates the satisfaction
        ! of the Kuhn-Tucker conditions
        if (wmax <= 0.0_PREC) exit

        iz = izmax
        j = index(iz)

        ! Step 5
        ! The sign of W(j) > 0 is OK for J to be moved to set P.
        ! Begin the transformation and check new diagonal element to avoid
        ! near linear dependence.

        tmp = a(np+1, j)
        call h12 (1, np+1, np+1+1, a(:,j), up)

        unorm = NRM2 (np, a(1:np,j), 1)

        diff = (unorm + abs(a(np+1,j)*FACTOR)) - unorm
        if (diff > 0.0_PREC) then
            ! Column J is sufficiently independent. Copy B into Z,
            ! update Z and solve for ZTEST (= proposed new value for X(j))
            z(1:m) = b(1:m)
            call h12 (2, np+1, np+1+1, a(1:m,j), up, z(1:m))

            ztest = z(np+1) / a(np+1, j)
        end if

        if (diff <= 0.0_PREC .or. ztest <= 0.0_PREC) then
            ! Reject J as a candidate to be moved from set Z to set P.
            ! Restore A(NPP1,J), set W(J) = 0 and loop back to test
            ! dual again.
            ! Note that if diff <= 0 we never compute ZTEST, but also
            ! don't care about its value.
            a(np+1,j) = tmp
            w(j) = 0.0_PREC
            cycle main
        end if

        ! The index J = INDEX(IZ) has been selected to be moved from
        ! set Z to set P. Update B, update indices, apply Householder
        ! transformations to columns in new set Z, zero subdiagonal
        ! elements in column J, set w(j) = 0.
        b(1:m) = z(1:m)

        index(iz) = index(iz1)
        index(iz1) = j
        iz1 = iz1 + 1
        np = np+1

        ! Apply Householder transformations to columns in new set Z
        do jz = iz1, iz2
            jj = index(jz)
            call h12 (2, np, np+1, a(:,j), up, a(:,jj))
        end do

        ! Zero out columns in A_p which do not correspond to some J in P
        a(np+1:m,j) = 0.0_PREC
        w(j) = 0.0_PREC

        ! Step 6: Solve the triangular system, temporarily store solution in Z
        call solve_triangular (a, np, index, z)

        ! === Secondary loop ===

        ! This loops over steps 6-11

        ! Theoretical upper bound on required number of iterations is 3*N
        loop6: do iter = 1, 3 * n
            ! See if all new constrained coefficients are feasible

            ! Note: in the original F77 code JJ was set in the "internal"
            ! subroutine here replaced by SOLVE_TRIANGULAR.
            ! Set it to same value as it would have had after exiting that
            ! code block.
            jj = index(1)

            ! Step 8+9: Find index q (stored here in JJ) such that
            ! X(q)/(X(q) - Z(q)) = min{X(i)/(X(i) - Z(i)) | z(i) <= 0, i in P}

            ! Set ALPHA to a value that will not be attained if
            ! Z(j) > 0.0 for all j in P.
            alpha = 2.0_PREC

            do ip = 1, np
                i = index(ip)
                if (z(ip) <= 0.0_PREC) then
                    t = x(i) / (x(i) - z(ip))
                    ! Keep track of minimal T
                    if (t < alpha) then
                        alpha = t
                        jj = ip
                    end if
                end if
            end do

            ! Step 7: Check that Z(J) > 0 for all J in P.
            ! If all new constrained coefficients are feasible alpha will
            ! still be 2.0. If so, exit from secondary loop to main loop.
            if (alpha >= 2.0_PREC) exit loop6

            ! Step 10:
            ! Otherwise use ALPHA which will be in [0,1] to interpolate
            ! between the old X and the new Z.
            do ip = 1, np
                i = index(ip)
                x(i) = x(i) + alpha * (z(ip) - x(i))
            end do

            ! Step 11: Modify A and B and the INDEX array to move coefficient
            ! JJ from set P to set Z
            loop11: do while (np > 0)
                i = index(jj)
                x(i) = 0.0_PREC

                if (jj /= np) then
                    jj = jj + 1
                    call rotate (jj, np, index, a, b)
                end if

                np = np - 1
                iz1 = iz1 - 1
                index(iz1) = i

                ! See if the remaining coefficients in set P are feasible.
                ! They should be because of the way alpha was determined.
                ! If any are infeasible it is due to round-off error.
                ! Any that are nonpositive will be set to zero and moved
                ! from set P to set Z.

                do ip = 1, np
                    i = index(ip)
                    if (x(i) <= 0.0_PREC) cycle loop11
                end do

                ! Exit loop, no non-positive elements in X found
                if (ip > np) exit
            end do loop11

            ! Copy B into Z, solve triangular system again and cycle
            ! loop for algorithm steps 6-11
            z(1:m) = b(1:m)
            call solve_triangular (a, np, index, z)

        end do loop6

        ! Secondary loop exceeded iteration limit
        if (iter > 3*n) then
            lstatus = NF_STATUS_MAX_ITER
            ! Original status code used in Lawson/Hanson
            lstatus%code_orig = STATUS_MAX_ITER
            if (present(msg)) msg = 'NNLS: Max. iteration count exceeded'
            exit main
        end if

        ! Regular termination of secondary loop within max. iteration limit
        do ip = 1, np
            i = index(ip)
            x(i) = z(ip)
        end do

        ! All new coefficients are positive. Cycle main loop

        ! Repeat this computation here instead of at the beginning
        ! of the main loop, as the original F77 code cycles the main loop
        ! but skips this part.

        ! Step 2: Compute coefficients of the dual (negative gradient) vector W
        do iz = iz1, iz2
            j = index(iz)
            sm = sum(a(np+1:m,j) * b(np+1:m))
            w(j) = sm
        end do
    end do main

    if ((np + 1) <= m) then
        rnorm = NRM2 (m-np, b(np+1:m), 1)
    else
        rnorm = 0.0_PREC
        w(1:n) = 0.0_PREC
    end if

100 continue

    ! Write back ISTATUS flag
    if (present(status)) status = lstatus

    contains

    subroutine assert_input (cond, lmsg)
        logical, intent(in) :: cond
        character (*), intent(in) :: lmsg

        lstatus = NF_STATUS_OK
        lstatus%code_orig = STATUS_OK

        if (.not. cond) then
            lstatus = NF_STATUS_INVALID_ARG
            ! Original status code for invalid dimensions
            lstatus%code_orig = STATUS_INVALID_DIMS

            if (present(msg)) then
                msg = 'NNLS: Invalid input: ' // lmsg
            end if
        end if
    end subroutine

    pure subroutine solve_triangular (a, np, index, z)
        real (PREC), intent(in), dimension(:,:), contiguous :: a
        integer, intent(in) :: np
        integer, intent(in), dimension(:), contiguous :: index
        real (PREC), intent(inout), dimension(:), contiguous :: z

        integer :: i, j, ip

        j = 0

        do ip = np, 1, -1
            if (ip < np) then
                ! Skip this part in the first iteration
                do i = 1, ip
                    z(i) = z(i) - a(i,j) * z(ip + 1)
                end do
            end if

            j = index(ip)
            z(ip) = z(ip) / a(ip, j)
        end do

    end subroutine


    subroutine rotate (jj, np, index, a, b)
        integer, intent(in) :: jj, np
        integer, intent(inout), dimension(:), contiguous :: index
        real (PREC), intent(inout), dimension(:,:), contiguous :: a
        real (PREC), intent(inout), dimension(:), contiguous :: b

        integer :: ip, l, i, n
        real (PREC) :: tmp, c, s, r

        n = size(a, 2)

        do ip = jj, np
            i = index(ip)
            index(ip-1) = i

            ! Call LAPACK routine to compute parameters of
            ! Givens rotation. This replaces the call to G1()
            ! in Lawson/Hanson. Note that we need to use LARTGP
            ! and not LARTG, as the latter flips the signs in
            ! c, s and r.
            call LARTGP (a(ip-1,i), a(ip,i), c, s, r)
            a(ip-1,i) = r
            a(ip,i) = 0.0_PREC

            do l = 1, n
                if (l /= i) then
                    ! Apply rotation to elements in A
                    tmp =  a(ip-1,l)
                    a(ip-1,l) =  c * tmp + s * a(ip,l)
                    a(ip,l)   = -s * tmp + c * a(ip,l)
                end if
            end do

            ! Apply rotation to elements in B
            tmp = b(ip-1)
            b(ip-1) =  c * tmp + s * b(ip)
            b(ip)   = -s * tmp + c * b(ip)
        end do

    end subroutine rotate


end subroutine



subroutine ldp_query (m, n, w, index)
    integer, intent(in) :: m, n
    integer, intent(out) :: w, index

    integer :: nnls_w, nnls_z, nnls_index
    ! Determine workspace requirements by NNLS
    call nnls_query (n+1, m, nnls_w, nnls_z, nnls_index)

    ! Total workspace requirement
    ! The working space W needs to store:
    !   1.  (N+1)*M elements = Matrix A for NNLS
    !   2.  N+1     elements = vector B for NNLS
    !   3.  M       elements = vector X for NNLS
    !   4.  Workspace required for W and Z arguments of NNLS
    w = (n+1) * (m+1) + m + nnls_w + nnls_z

    index = nnls_index

end subroutine



subroutine ldp (g, h, x, xnorm, w, index, status, msg)
    !*  LDP solves the constrainted least-distance problem
    !       min.    ||x||   s.t. Gx >= h
    !   using Algorithm 23.27 in
    !       Lawson, Hanson (1995): Solving least squares problems
    !
    !   Call subroutine NNLS to do the actual work.
    real (PREC), intent(inout), dimension(:,:), contiguous :: g
        !*  Matrix G with dimension (M,N).
    real (PREC), intent(in), dimension(:), contiguous :: h
        !*  Array H contains the vector h of the problem.
    real (PREC), intent(out), dimension(:), contiguous :: x
        !*  Contains the vector X on normal termination (IMODE = 1).
        !   On an error termination, contents are set to zero.
    real (PREC), intent(out) :: xnorm
        !*  On normal termination, XNORM contains the Euclidean norm of X.
        !   On error terminal, XNORM is set to zero.
    real (PREC), intent(out), dimension(:), contiguous, target :: w
        !*  Working space with dimension of at least (M + 2)*(N + 1) + 2*M
    integer, intent(out), dimension(:), contiguous :: index
        !*  Working space with dimension of at least M.
    type (status_t), intent(out), optional :: status
        !*  Optional exit status flag. If present and
        !   ISTATUS = STATUS_WORKSPACE_QUERY on input, routine writes
        !   required workspace array sizes into the first elements of
        !   W and IP and immediately exits.
    character (*), intent(out), optional :: msg
        !*  Optional status message. Changed by routine only on
        !   unsuccessful termination.

    integer :: m, n, j
    real (PREC) :: fac

    ! Variables used to call NNLS
    integer :: ia, ib, ix, iz, iw
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_a
    real (PREC), dimension(:), pointer, contiguous :: ptr_b, ptr_z, ptr_x, ptr_w
    real (PREC) :: rnorm
    integer :: nnls_nw, nnls_nz, nnls_nindex, nw, nindex

    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK
    ! Original success exit code
    lstatus%code_orig = STATUS_OK

    m = size(g, 1)
    n = size(g, 2)

    ! === Workspace dimensions ===

    ! Determine workspace requirements by NNLS
    call nnls_query (n+1, m, nnls_nw, nnls_nz, nnls_nindex)

    ! Total workspace requirement
    call ldp_query (m, n, nw, nindex)

    ! === Input checks ===

    call assert_input (size(g) > 0, 'Array G has invalid size')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(h) == m, 'Non-conformable arrays G, H')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(x) == n, 'Non-conformable arrays G, X')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(w) >= nw, 'Workspace W too small')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(index) >= nindex, 'INDEX array too small')
    if (lstatus /= NF_STATUS_OK) goto 100

    ! === LDP implementation ===

    xnorm = 0.0_PREC
    x = 0.0_PREC

    ! === Copy data into arrays used to call NNLS ===

    ! Indices on working space W
    ia = 1
    ib = ia + (n+1)*m
    iz = ib + (n+1)
    ix = iz + nnls_nz
    iw = ix + m

    ! Pointers to blocks of W
    ptr_A(1:n+1,1:m) => w(1:(n+1)*m)
    ptr_b(1:n+1) => w(ib:ib+n)
    ptr_z(1:nnls_nz) => w(iz:iz+nnls_nz-1)
    ptr_x(1:m) => w(ix:ix+m-1)
    ptr_w(1:nnls_nw) => w(iw:iw+nnls_nw-1)

    ! Copy G^T into A
    do j = 1, m
        ptr_A(1:n,j) = g(j,:)
    end do
    ! Copy h^T
    ptr_A(n+1,:) = h

    ! vector b of Ax = b for NNLS
    ptr_b(:) = 0.0
    ptr_b(n+1) = 1.0

    ! === LPD solver ===

    ! Compute solution to dual problem:
    ! Use NNLS to compute augmented problem
    !   Ax = b  s.t. x >= 0
    ! where A^T = [G^T, h^T], b^T = [0,...,0,1]^T
    call nnls (ptr_a, ptr_b, ptr_x, rnorm, ptr_w, ptr_z, index, lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (rnorm <= 0.0_PREC) then
        ! Incompatible constraints
        lstatus = NF_STATUS_INVALID_STATE
        lstatus%code_orig = STATUS_INCOMPAT_CONSTR
        if (present(msg)) msg = 'LDP: Incompatible inequality constraints'
        goto 100
    end if

    ! Step 2: Compute residual vector
    fac = 1.0_PREC - dot_product (ptr_x(1:m), h(1:m))

    if (abs(fac) < epsilon(1.0_PREC)) then
        lstatus = NF_STATUS_INVALID_STATE
        lstatus%code_orig = STATUS_INCOMPAT_CONSTR
        if (present(msg)) msg = 'LDP: Incompatible inequality constraints'
        goto 100
    end if

    ! Compute solution to primal problem:
    ! Use solution vector from NNLS compute solution vector to LDP
    fac = 1.0_PREC / fac
    do j = 1, n
        x(j) = dot_product (g(1:m,j), ptr_x(1:m))
    end do

    x = x * fac

    xnorm = NRM2 (n, x, 1)

    ! Compute Lagrange multipliers for primal problem
    ! (this part is not present in original Lawson/Hanson code, but
    ! comes from Dieter Kraft's implementation -- required for SLSQP).
    do j = 1, m
        w(j) = fac * ptr_x(j)
    end do

100 continue

    if (present(status)) status = lstatus

    contains

    subroutine assert_input (cond, lmsg)
        logical, intent(in) :: cond
        character (*), intent(in) :: lmsg

        lstatus = NF_STATUS_OK
        lstatus%code_orig = STATUS_OK
        if (.not. cond) then
            lstatus = NF_STATUS_INVALID_ARG
            ! Original status code for invalid dimensions
            lstatus%code_orig = STATUS_INVALID_DIMS

            if (present(msg)) then
                msg = 'LDP: Invalid input: ' // lmsg
            end if
        end if
    end subroutine

end subroutine
