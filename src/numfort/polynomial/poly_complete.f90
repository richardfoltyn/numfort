

module numfort_polynomial_complete
    !*  Module contains some utility functions for complete polynomials
    !   that are common to all kind/type-specific routines handling
    !   complete polynomials.

    use, intrinsic :: iso_fortran_env

    use numfort_core
    use numfort_common
    use numfort_common_input_checks

    implicit none
    private

    public :: poly_complete_get_nterms
    public :: polyexponents_complete

    contains




elemental function poly_complete_get_nterms (n, k) result(res)
    !*  POLY_COMPLETE_GET_NTERMS returns the number of terms in a complete
    !   polynomial in N variables of degree K.
    integer, intent(in) :: n
        !*  Number of variables
    integer, intent(in) :: k
        !*  (Maximum) degree of complete polynomial
    integer :: res
        !*  Number of terms

    ! Number of terms is obtained as follows: each basis function is drawn from
    ! the set (1,x1,...xn) as the product of there elements.
    ! Elements can be drawn multiple times, but order does not matter as the
    ! product is unchanged.
    ! Hence this is an application of choosing k elements from a set of n
    ! elements with repetition, but ignoring the order, and the number of
    ! terms is given by binom{n+k-1}{k}
    res = comb (n+1, k, repetition=.true.)

end function



pure subroutine polyexponents_complete (n, k, exponents, status)
    !*  POLYEXPONENTS_COMPLETE returns the list of all exponent permutations
    !   that occur in a complete polynomial of given degree for a given
    !   number of variables.
    !
    !   Note on the implementation:
    !
    !   The exponent permutations are computed blockwise, such that for each
    !   "block" the exponents sum up to the same fixed number k_i in [0,..,k].
    !   The first two blocks are thus given by
    !
    !       [0,0,.....,0]
    !
    !       [1,0,......0]
    !       [0,1,0,...,0]
    !       [...........]
    !       [0,.....,0,1]
    !
    !   ie they are the zero n-vector and the n-by-n identity matrix.
    !   Any consecutive blocks are formed from the previous block by
    !       1.  Copying the entire previous block and adding 1 to all exponents
    !           for variable 1, x_1.
    !       2.  Copying those rows from the previous block where the
    !           exponent for x_1 is zero, and adding 1 to all exponents for
    !           x_2.
    !       3.  etc. until we process x_n
    !   We need to restrict which rows are copied and incremented for all
    !   variables other than x_1 to avoid introducing duplicates.

    integer, intent(in) :: n
        !*  Number of variables
    integer, intent(in) :: k
        !*  (Maximum) degree of complete polynomial
    integer, intent(out), dimension(:,:), contiguous :: exponents
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: nterms, ioffset, ioffset_prev, ni, i, j, ic
    integer, dimension(:), allocatable :: nterms_i

    lstatus = NF_STATUS_OK

    ! Input checks
    call check_nonneg (n, n, "n", lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100
    call check_nonneg (k, k, "k", lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! compute total number of terms
    nterms = poly_complete_get_nterms (n, k)

    if (size(exponents,2) < nterms) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Constant term
    exponents(:,1) = 0

    if (k == 0 .or. n == 0) goto 100

    ! For each "block" with maximum degree j in {1,..,k}, keep track of
    ! how many of the columns in PREVIOUS block were created by increasing
    ! the exponent of variable i in {1,...,n}
    allocate (nterms_i(n), source=0)

    ioffset_prev = 0
    ioffset = 1

    ! Loop over remaining polynomial degrees j in {1,...,k}
    do j = 1, k

        ! loop over each variable
        do i = 1, n
            ! Take max to ensure that there is at least one term to be
            ! processed for variable i, which corresponds to the monomial
            ! x_i^j and needs to be obtained from x_i^(j-1) from the previous
            ! block.
            ! Note: this is only an issue for j=1 and not necessary for
            ! higher degrees.
            ni = max(sum(nterms_i(i:n)), 1)
            ! Offset within the previous block corresponding to degree (j-1)
            ! When processing variable i, we ignore all terms in columns
            ! prior to IOFFSET_PREV as this would create duplicates.
            if (i > 1) then
                ioffset_prev = ioffset_prev + nterms_i(i-1)
            end if
            do ic = 1, ni
                exponents(:,ioffset+ic) = exponents(:,ioffset_prev+ic)
                ! Increment exponent on variable i
                exponents(i,ioffset+ic) = exponents(i,ioffset+ic) + 1
            end do
            ioffset = ioffset + ni
        end do

        ! Update stuff for next iteration

        ! Update NTERMS for each variable
        do i = 1, n
            ! Again, ensure that the min. number of terms is 1.
            nterms_i(i) = max(sum(nterms_i(i:n)), 1)
        end do

        ! Reset prev. offset to the *beginning* of block that we just processed
        ioffset_prev = comb (n+1, j-1, repetition=.true.)
        ! Set offset to beginning of next block
        ioffset = ioffset_prev + sum(nterms_i)
    end do

100 continue

    if (allocated(nterms_i)) deallocate (nterms_i)

    if (present(status)) status = lstatus

end subroutine

end module