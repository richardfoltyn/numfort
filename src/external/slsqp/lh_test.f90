

program lh_test

    use, intrinsic :: iso_fortran_env, only: real64

    use lawson_hanson_orig_real64, h12_orig => h12, nnls_orig => nnls, &
        ldp_orig => ldp
    use lawson_hanson_real64

    integer, parameter :: PREC = real64

    call test_h12 ()
    call test_nnls ()
    call test_ldp ()

    contains


subroutine test_h12 ()

    integer, parameter :: n = 5
    integer :: i, j

    real (PREC), dimension(n,n) :: a, aorig
    real (PREC), dimension(n) :: up
    real (PREC) :: diff


    do j = 1, n
        do i = 1, n
            a(i,j) = (j-1) * n + i
        end do
    end do

    aorig = a

    do j = 3, n
        call h12_orig (1,j-1,j,n,aorig(1,j-2),1,up(j-2),aorig(1,j-1),1,n,n-j+2)
        call h12 (1, j-1, j, a(:,j-2), up(j-2), a(:,j-1:n))

        diff = maxval(abs(a-aorig))

        call h12_orig (2,j-1,j,n,aorig(1,j-2),1,up(j-2),aorig,n,1,n)
        call h12 (2, j-1, j, a(:,j-2), up(j-2), a, trans=.true.)

        diff = maxval(abs(a-aorig))

    end do

end subroutine



subroutine test_nnls ()

    integer, parameter :: m = 10
    integer, parameter :: n = 5

    real (PREC), dimension(m,n) :: a, aorig
    real (PREC), dimension(m) :: b, borig
    real (PREC), dimension(max(m,n)) :: w, worig
    real (PREC), dimension(n) :: x, xorig
    real (PREC), dimension(max(m,n)) :: z, zorig
    integer, dimension(n) :: index, index_orig

    real (PREC) :: rnorm, rnorm_orig
    integer :: mode, mode_orig, istatus
    character (100) :: msg

    integer :: rnd_size
    integer, allocatable, dimension(:) :: rnd_put

    iprint = 1

    call random_seed (size=rnd_size)
    allocate (rnd_put(rnd_size))
    rnd_put(:) = 1
    call random_seed (put=rnd_put)

    call random_number (a)

    a(:,1) = 1.0
    a(:,2) = a(:,2) * 10.0
    a(:,3) = a(:,2)
    a(:,5) = a(:,2) + a(:,4)

    x = [3.0d0, 0.5d0, 1.0d0, 2.0d0, 0.5d0]
    call dgemv ('N', m, n, 1.0d0, a, m, x, 1, 0.0d0, b, 1)

    x = 0.0

    aorig = a
    borig = b

    x = 0.0
    w = 0.0
    z = 0.0
    index = 0

    xorig = 0.0
    worig = 0.0
    zorig = 0.0
    index_orig = 0


    call nnls_orig (aorig, m, m, n, borig, xorig, rnorm_orig, worig, zorig, index_orig, mode_orig)
    call nnls (a, b, x, rnorm, w, z, index, istatus, msg)

    print *, x

    print 100, 'A', maxval(abs(a-aorig))
    print 100, 'b', maxval(abs(b-borig))
    print 100, 'w', maxval(abs(w-worig))
    print 100, 'z', maxval(abs(z-zorig))
    print 100, 'index', real(maxval(abs(index-index_orig)), real64)
    print 100, 'rnorm', rnorm-rnorm_orig

100 format ('Diff ', /a, ': ', es10.2e2)

end subroutine



subroutine test_ldp ()
    !*  Test routine for the minimum distance problem
    !       min ||x||   s.t. Gx >= h
    integer, parameter :: m = 10
    integer, parameter :: n = 5

    real (PREC), dimension(m,n) :: g, gorig
    real (PREC), dimension(m) :: h, horig
    real (PREC), dimension((m+2)*(n+1)+2*m) :: worig
    real (PREC), dimension(:), allocatable :: w
    real (PREC), dimension(n) :: x, xorig
    integer, dimension(m) :: index, index_orig
    integer :: nw, nindex

    real (PREC) :: xnorm, xnorm_orig
    integer :: mode_orig, istatus
    character (100) :: msg

    integer :: rnd_size
    integer, allocatable, dimension(:) :: rnd_put

    iprint = 1

    call random_seed (size=rnd_size)
    allocate (rnd_put(rnd_size))
    rnd_put(:) = 2
    call random_seed (put=rnd_put)

    call random_number (g)

    x = [3.0d0, 0.5d0, 1.0d0, 2.0d0, 0.5d0]
    call random_number (h)

    call dgemv ('N', m, n, 1.0d0, g, m, x, 1, -1.0d0, h, 1)

    x = 0.0

    gorig = g
    horig = h

    x = 0.0
    index = 0

    xorig = 0.0
    worig = 0.0
    index_orig = 0


    call ldp_orig (gorig, m, m, n, horig, xorig, xnorm_orig, worig, index_orig, mode_orig)

    call ldp_query (size(g, 1), size(g, 2), nw, nindex)
    allocate (w(nw), source=0.0_PREC)
    call ldp (g, h, x, xnorm, w, index, istatus, msg)

    print *, x

    print 100, 'A', maxval(abs(g-gorig))
    print 100, 'x', maxval(abs(x-xorig))
    print 100, 'w', maxval(abs(w-worig))
    print 100, 'index', real(maxval(abs(index-index_orig)), real64)
    print 100, 'xnorm', xnorm-xnorm_orig

100 format ('Diff ', /a, ': ', es10.2e2)


end subroutine




end program
