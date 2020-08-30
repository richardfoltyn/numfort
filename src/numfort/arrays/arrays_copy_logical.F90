

module numfort_arrays_copy_logical

    use numfort_common_status
    use numfort_common_input_checks

    implicit none

    private

    public :: pack_indexed

    interface pack_indexed
        procedure pack_indexed_dim_2d
    end interface


    contains



subroutine pack_indexed_dim_2d (src, indices, dst, dim, trans, status)
    !*  PACK_INDEXED copies a subset of elements from SRC to DST, where
    !   an elements of DST is given by
    !       DST(i,j) = SRC(INDICES(i),j)    if DIM=1 and TRANS=.FALSE.
    !       DST(i,j) = SRC(INDICES(j),i)    if DIM=1 and TRANS=.TRUE.
    !       DST(i,j) = SRC(i,INDICES(j))    if DIM=2 and TRANS=.FALSE.
    !       DST(i,j) = SRC(j,INDICES(i))    if DIM=2 and TRANS=.TRUE.
    logical, intent(in), dimension(:,:), contiguous :: src
    integer, intent(in), dimension(:), contiguous :: indices
    logical, intent(out), dimension(:,:), contiguous :: dst
    integer, intent(in), optional :: dim
    logical, intent(in), optional :: trans
    type (status_t), intent(out), optional :: status

    integer :: ldim
    logical :: ltrans
    integer :: i, k, dim_dst
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'COPY_INDEXED'

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (dim, 1, ldim)
    call set_optional_arg (trans, .false., ltrans)

    call check_cond (ldim == 1 .or. ldim == 2, NAME, 'DIM: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if ((ltrans .and. ldim == 1) .or. (.not. ltrans .and. ldim == 2)) then
        dim_dst = 2
    else
        dim_dst = 1
    end if

    call check_cond (size(dst,dim_dst) == size(indices) .and. &
        size(dst,3-dim_dst) == size(src,3-ldim), NAME, &
        'DST: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Copy data ---

    if (ldim == 1) then
        if (ltrans) then
            do i = 1, size(indices)
                k = indices(i)
                dst(:,i) = src(k,:)
            end do
        else
            do i = 1, size(indices)
                k = indices(i)
                dst(i,:) = src(k,:)
            end do
        end if
    else
        if (ltrans) then
            do i = 1, size(indices)
                k = indices(i)
                dst(i,:) = src(:,k)
            end do
        else
            do i = 1, size(indices)
                k = indices(i)
                dst(:,i) = src(:,k)
            end do
        end if
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine


end module
