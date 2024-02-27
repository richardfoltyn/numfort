

subroutine __APPEND(test_gemm,__PREC) (tests)
    !* TEST_GEMM implements unit tests for the F95 wrapper for GEMM.
    use blas_interfaces, only: BLAS_GEMM => GEMM

    integer, parameter  :: PREC = __PREC
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (str) :: msg

    integer :: m, n, k, lda, ldb, ldc
    integer, parameter :: mrange(*) = [1, 2, 10]
    integer, parameter :: nrange(*) = [1, 2, 10]
    integer, parameter :: krange(*) = [1, 2, 10]
    real (PREC), dimension(:,:), allocatable :: a, b, c, c95
    real (PREC) :: alpha, beta
    character (1) :: transa, transb
    character (1) :: trans_range(2) = ['N', 'T']

    integer :: i, j, ik, ita, itb

    tc => tests%add_test ("GEMM wrapper for real*" // str(PREC))

    call set_seed (123)

    do ita = 1, size(trans_range)
        transa = trans_range(ita)
        do itb = 1, size(trans_range)
            transb = trans_range(itb)

            do i = 1, size(mrange)
                m = mrange(i)

                do j = 1, size(nrange)
                    n = nrange(j)

                    do ik = 1, size(krange)
                        k = krange(ik)

                        if (transa == 'N') then
                            allocate (a(m, k))
                            lda = m
                        else
                            allocate (a(k, m))
                            lda = k
                        end if

                        if (transb == 'N') then
                            allocate (b(k, n))
                            ldb = k
                        else
                            allocate (b(n, k))
                            ldb = n
                        end if

                        allocate (c(m, n), c95(m, n))
                        ldc = m

                        call random_number (a)
                        call random_number (b)
                        call random_number (c)

                        c95(:,:) = c

                        ! === Call with default ALPHA, BETA ===
                        alpha = 1.0
                        beta = 0.0

                        ! Call F77 interface
                        call BLAS_GEMM (transa, transb, m, n, k, alpha, a, lda, &
                            b, ldb, beta, c, ldc)

                        ! Call F95 interface (with default alpha, beta) and default
                        ! TRANS{A,B} (if applicable)
                        if (transa == 'N' .and. transb == 'N') then
                            call GEMM (a, b, c95)
                        else
                            call GEMM (a, b, c95, transa=transa, transb=transb)
                        end if

                        msg = 'GEMM (transa=' // transa // 'transb=' // transb &
                            // ', m=' // str(m) // ', n=' // str(n) &
                            // ', k=' // str(k) // ')'
                        call tc%assert_true (all(c==c95), msg)

                        call GEMM (a, b, c95, alpha=alpha, beta=beta, &
                            transa=transa, transb=transb)

                        msg = 'GEMM (transa=' // transa // 'transb=' // transb &
                            // ', m=' // str(m) // ', n=' // str(n) &
                            // ', k=' // str(k) &
                            // ', alpha=' // str(alpha, 'f0.2') // ')'
                        call tc%assert_true (all(c==c95), msg)

                        ! === Call with non-default ALPHA, BETA ===
                        c95(:,:) = c
                        call random_number (alpha)
                        call random_number (beta)
                        ! Call F77 interface
                        call BLAS_GEMM (transa, transb, m, n, k, alpha, a, lda, &
                            b, ldb, beta, c, ldc)
                        call GEMM (a, b, c95, alpha=alpha, beta=beta, &
                            transa=transa, transb=transb)

                        msg = 'GEMM (transa=' // transa // 'transb=' // transb &
                            // ', m=' // str(m) // ', n=' // str(n) &
                            // ', k=' // str(k) &
                            // ', alpha=' // str(alpha, 'f0.2') // ')'
                        call tc%assert_true (all(c==c95), msg)

                        deallocate (a, b, c, c95)

                    end do
                end do
            end do
        end do
    end do

end subroutine
