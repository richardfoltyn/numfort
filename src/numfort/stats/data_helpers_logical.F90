

module numfort_stats_data_helpers_logical

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays, only: arange, unique
    use numfort_common_alloc
    use numfort_common_cond_alloc
    use numfort_common_status
    use numfort_common_input_checks
    use numfort_stats_core

    implicit none
    private

    public :: get_regressor_dims
    public :: extract_block
    public :: extract_block_alloc

    interface get_regressor_dims
        procedure get_regressor_dims
    end interface

    interface extract_block
        procedure extract_block_1d, extract_block_2d
    end interface

    interface extract_block_alloc
        procedure extract_block_alloc_1d, extract_block_alloc_2d
    end interface

    contains


pure subroutine get_regressor_dims (X, trans, nvars, nobs)
    !*  GET_REGRESSOR_DIMS returns the dimensions of the matrix of regressors,
    !   depending on whether the matrix is provided in transposed form or not.
    logical, intent(in), dimension(:,:) :: X
    logical, intent(in), optional :: trans
    !*  If true, regressor matrix is provided in transposed form
    !   (variables along dimension 1, observations along dimension 2)
    integer, intent(out), optional :: nvars
    !*  Number of variables
    integer, intent(out), optional :: nobs
    !*  Number of observations

    logical :: ltrans
    integer :: lnobs, lnvars, dim

    ltrans = .false.
    if (present(trans)) ltrans = trans

    ! Dimension along with observations are aligned
    dim = 1
    if (ltrans) dim = 2

    lnobs = size(X, dim)
    lnvars = size(X, 3-dim)

    if (present(nvars)) nvars = lnvars
    if (present(nobs)) nobs = lnobs

end subroutine


subroutine extract_block_2d (x, ifrom, n, trans, x_block, x_rest, status)
    !*  EXTRACT_BLOCK splits a contiguous block from a given array X,
    !   and optionally returns this block and/or the remainder of the array.
    logical, intent(in), dimension(:,:), contiguous :: x
        !*  Array to be split
    integer, intent(in) :: ifrom
        !*  Initial index of the block to be extracted. The index is assumed
        !   to correspond to the first dimension, unless TRANS=.TRUE.
        !   is specified.
    integer, intent(in) :: n
        !*  Size of the contiguous block to be extracted.
    logical, intent(in), optional :: trans
        !*  If present and true, store data in transposed form in output
        !   arrays (default: false)
    logical, intent(out), dimension(:,:), contiguous, optional :: x_block
        !*  Array containing the contiguous block extracted from X.
    logical, intent(out), dimension(:,:), contiguous, optional :: x_rest
        !*  Array containing the remained of X, with X_BLOCK removed.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: i, Nvars, Nobs, Nrest
    logical :: ltrans
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'EXTRACT_BLOCK'

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (trans, .false., ltrans)

    call get_regressor_dims (x, ltrans, Nvars, Nobs)
    Nrest = Nobs - n

    if (present(x_block)) then
        call check_cond (size(x_block,1) == n .and. size(x_block,2) == Nvars, &
            NAME, 'X_BLOCK: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(x_rest)) then
        call check_cond (size(x_rest,1) == Nrest .and. size(x_rest,2) == Nvars, &
            NAME, 'X_REST: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Copy data ---

    if (present(x_block)) then
        if (ltrans) then
            do i = 1, Nvars
                x_block(1:n,i) = x(i,ifrom:ifrom+n-1)
            end do
        else
            do i = 1, Nvars
                x_block(1:n,i) = x(ifrom:ifrom+n-1,i)
            end do
        end if
    end if

    if (present(x_rest)) then
        if (ltrans) then
            do i = 1, Nvars
                x_rest(1:ifrom-1,i) = x(i,1:ifrom-1)
                x_rest(ifrom:Nrest,i) = x(i,ifrom+n:Nobs)
            end do
        else
            do i = 1, Nvars
                x_rest(1:ifrom-1,i) =  x(1:ifrom-1,i)
                x_rest(ifrom:Nrest,i) = x(ifrom+n:Nobs,i)
            end do
        end if
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine extract_block_alloc_2d (x, ifrom, n, trans, x_block, x_rest, status)
    !*  EXTRACT_BLOCK_ALLOC splits a contiguous block from a given array X,
    !   and optionally returns this block and/or the remainder of the array.
    !
    !   The arrays X_BLOCK and X_REST are (re)allocated to store the
    !   data if necessary.
    logical, intent(in), dimension(:,:), contiguous :: x
        !*  Array to be split
    integer, intent(in) :: ifrom
        !*  Initial index of the block to be extracted. The index is assumed
        !   to correspond to the first dimension, unless TRANS=.TRUE.
        !   is specified.
    integer, intent(in) :: n
        !*  Size of the contiguous block to be extracted.
    logical, intent(in), optional :: trans
        !*  If present and true, store data in transposed form in output
        !   arrays (default: false)
    logical, intent(inout), dimension(:,:), allocatable, optional :: x_block
        !*  Array containing the contiguous block extracted from X.
        !   If the input array does not have the exact size required to
        !   hold the data, it is (re)allocated.
    logical, intent(inout), dimension(:,:), allocatable, optional :: x_rest
        !*  Array containing the contiguous block extracted from X.
        !   If the input array does not have the exact size required to
        !   hold the data, it is (re)allocated.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nvars, Nobs, Nrest
    integer :: shp(2)

    call get_regressor_dims (x, trans, Nvars, Nobs)
    Nrest = Nobs - n

    if (present(x_block)) then
        shp(1) = n
        shp(2) = Nvars
        call cond_alloc (x_block, shp)
    end if

    if (present(x_rest)) then
        shp(1) = Nrest
        shp(2) = Nvars
        call cond_alloc (x_rest, shp)
    end if

    ! ifort-2019 seems to have troubles handling optional, allocatable
    ! dummy arguments when passing them to routines. Split out calls
    ! depending on which arguments are actually present.
    if (present(x_block) .and. present(x_rest)) then
        call extract_block (x, ifrom, n, trans, x_block, x_rest, status)
    else if (present(x_block)) then
        call extract_block (x, ifrom, n, trans, x_block=x_block, status=status)
    else if (present(x_rest)) then
        call extract_block (x, ifrom, n, trans, x_rest=x_rest, status=status)
    end if

end subroutine



subroutine extract_block_1d (x, ifrom, n, x_block, x_rest, status)
    !*  EXTRACT_BLOCK splits a contiguous block from a given array X,
    !   and optionally returns this block and/or the remainder of the array.
    logical, intent(in), dimension(:), contiguous :: x
        !*  Array to be split
    integer, intent(in) :: ifrom
        !*  Initial index of the block to be extracted.
    integer, intent(in) :: n
        !*  Size of the contiguous block to be extracted.
    logical, intent(out), dimension(:), contiguous, optional :: x_block
        !*  Array containing the contiguous block extracted from X.
    logical, intent(out), dimension(:), contiguous, optional :: x_rest
        !*  Array containing the remained of X, with X_BLOCK removed.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nobs, Nrest
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'EXTRACT_BLOCK'

    lstatus = NF_STATUS_OK

    Nobs = size(x)
    Nrest = Nobs - n

    if (present(x_block)) then
        call check_cond (size(x_block) == n, NAME, &
            'X_BLOCK: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(x_rest)) then
        call check_cond (size(x_rest) == Nrest, NAME, &
            'X_REST: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Copy data ---

    if (present(x_block)) then
        x_block(1:n) = x(ifrom:ifrom+n-1)
    end if

    if (present(x_rest)) then
        x_rest(1:ifrom-1) =  x(1:ifrom-1)
        x_rest(ifrom:Nrest) = x(ifrom+n:Nobs)
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine extract_block_alloc_1d (x, ifrom, n, x_block, x_rest, status)
    !*  EXTRACT_BLOCK_ALLOC splits a contiguous block from a given array X,
    !   and optionally returns this block and/or the remainder of the array.
    !
    !   The arrays X_BLOCK and X_REST are (re)allocated to store the
    !   data if necessary.
    logical, intent(in), dimension(:), contiguous :: x
        !*  Array to be split
    integer, intent(in) :: ifrom
        !*  Initial index of the block to be extracted.
    integer, intent(in) :: n
        !*  Size of the contiguous block to be extracted.
    logical, intent(inout), dimension(:), allocatable, optional :: x_block
        !*  Array containing the contiguous block extracted from X.
        !   If the input array does not have the exact size required to
        !   hold the data, it is (re)allocated.
    logical, intent(inout), dimension(:), allocatable, optional :: x_rest
        !*  Array containing the contiguous block extracted from X.
        !   If the input array does not have the exact size required to
        !   hold the data, it is (re)allocated.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nobs, Nrest

    Nobs = size(x)
    Nrest = Nobs - n

    if (present(x_block)) then
        call cond_alloc (x_block, n)
    end if

    if (present(x_rest)) then
        call cond_alloc (x_rest, Nrest)
    end if

    ! ifort-2019 seems to have troubles handling optional, allocatable
    ! dummy arguments when passing them to routines. Split out calls
    ! depending on which arguments are actually present.
    if (present(x_block) .and. present(x_rest)) then
        call extract_block (x, ifrom, n, x_block, x_rest, status)
    else if (present(x_block)) then
        call extract_block (x, ifrom, n, x_block=x_block, status=status)
    else if (present(x_rest)) then
        call extract_block (x, ifrom, n, x_rest=x_rest, status=status)
    end if

end subroutine



end module
