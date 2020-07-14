


subroutine pack_indexed_2d (src, irow, icol, dst, trans)
    !*  PACK_INDEXED copies a subset of elements from SRC to DST, where
    !   an elements of DST is given by
    !       DST(i,j) = SRC(IROW(i),ICOL(j))     if TRANS=.FALSE.
    !       DST(j,i) = SRC(IROW(i),ICOL(j))     if TRANS=.TRUE.
    real (PREC), intent(in), dimension(:,:), contiguous :: src
        !*  Array from which data should be copied.
    integer, intent(in), dimension(:), contiguous :: irow
        !*  Mapping of row indices.
    integer, intent(in), dimension(:), contiguous :: icol
        !*  Mapping of column indices.
    real (PREC), intent(out), dimension(:,:) :: dst
        !*  Array to store packed data.
    logical, intent(in), optional :: trans
        !*  If present and true, transpose the sub-array in SRC when copying
        !   to DST.

    logical :: ltrans
    integer :: i, j, ir, jc

    call set_optional_arg (trans, .false., ltrans)

    if (ltrans) then
        do j = 1, size(icol)
            jc = icol(j)
            do i = 1, size(irow)
                ir = irow(i)
                dst(j,i) = src(ir,jc)
            end do
        end do
    else
        do j = 1, size(icol)
            jc = icol(j)
            do i = 1, size(irow)
                ir = irow(i)
                dst(i,j) = src(ir,jc)
            end do
        end do
    end if
end subroutine



subroutine pack_indexed_dim_2d (src, indices, dst, dim, trans, status)
    !*  PACK_INDEXED copies a subset of elements from SRC to DST, where
    !   an elements of DST is given by
    !       DST(i,j) = SRC(INDICES(i),j)    if DIM=1 and TRANS=.FALSE.
    !       DST(i,j) = SRC(INDICES(j),i)    if DIM=1 and TRANS=.TRUE.
    !       DST(i,j) = SRC(i,INDICES(j))    if DIM=2 and TRANS=.FALSE.
    !       DST(i,j) = SRC(j,INDICES(i))    if DIM=2 and TRANS=.TRUE.
    real (PREC), intent(in), dimension(:,:), contiguous :: src
    integer, intent(in), dimension(:), contiguous :: indices
    real (PREC), intent(out), dimension(:,:), contiguous :: dst
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



subroutine pack_indexed_1d (src, indices, dst, status)
    !*  PACK_INDEXED copies a subset of elements in SRC to DST, where an
    !   element in DST is given by
    !       DST(i) = SRC(INDICES(i))
    real (PREC), intent(in), dimension(:), contiguous :: src
    integer, intent(in), dimension(:), contiguous :: indices
    real (PREC), intent(out), dimension(:), contiguous :: dst
    type (status_t), intent(out), optional :: status

    integer :: i
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'COPY_INDEXED'

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call check_cond (size(dst) == size(indices), NAME, &
        'DST: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    do i = 1, size(indices)
        dst(i) = src(indices(i))
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine unpack_indexed_2d (src, irow, icol, dst, trans, fill)
    !*  UNPACK_INDEXED unpacks an array into a (weakly) larger array
    !   using a given (rectangular) mapping of row and column indices.
    !   Any element in the destination array that is not present in the
    !   source array is assigned a constant value, if present.
    !
    !   The elements in DST are given by
    !       DST(IROW(i),ICOL(j)) = SRC(i,j)     if TRANS=.FALSE.
    !       DST(IROW(i),ICOL(j)) = SRC(j,i)     if TRANS=.TRUE.
    real (PREC), intent(in), dimension(:,:), contiguous :: src
    integer, intent(in), dimension(:), contiguous :: irow
        !*  Mapping of row indices.
    integer, intent(in), dimension(:), contiguous :: icol
        !*  Mapping of column indices.
    real (PREC), intent(out), dimension(:,:) :: dst
    logical, intent(in), optional :: trans
        !*  If present and true, transpose the sub-array in SRC when copying
        !   to DST.
    real (PREC), intent(in), optional :: fill
        !*  Fill value assigned to all elements in DST which are not present
        !   in SRC. If not present, missing elements will remain untouched.

    logical :: ltrans
    integer :: i, j, ir, jc

    call set_optional_arg (trans, .false., ltrans)

    if (present(fill)) then
        dst = fill
    end if

    if (ltrans) then
        do j = 1, size(icol)
            jc = icol(j)
            do i = 1, size(irow)
                ir = irow(i)
                dst(ir,jc) = src(j,i)
            end do
        end do
    else
        do j = 1, size(icol)
            jc = icol(j)
            do i = 1, size(irow)
                ir = irow(i)
                dst(ir,jc) = src(i,j)
            end do
        end do
    end if
end subroutine



subroutine unpack_indexed_dim_2d (src, indices, dim, dst, trans, fill)
    !*  UNPACK_INDEXED unpacks an array into a (weakly) larger array
    !   using a given (rectangular) mapping indices along a given dimension.
    !   Any element in the destination array that is not present in the
    !   source array is assigned a constant value, if present.
    real (PREC), intent(in), dimension(:,:), contiguous :: src
    integer, intent(in), dimension(:), contiguous :: indices
        !*  Mapping of indices.
    integer, intent(in) :: dim
        !*  Dimension along which to operate. If DIM = 1, then rows
        !   will be unpacked as DST(INDICES(i),:) = SRC(i,:), and if
        !   DIM = 2, then DST(:,INDICES(i)) = SRC(:,i).
    real (PREC), intent(out), dimension(:,:) :: dst
    logical, intent(in), optional :: trans
    real (PREC), intent(in), optional :: fill
        !*  Fill value assigned to all elements in DST which are not present
        !   in SRC. If not present, missing elements will remain untouched.

    integer :: i, j
    logical :: ltrans

    call set_optional_arg (trans, .false., ltrans)

    if (present(fill)) then
        dst = fill
    end if

    if (dim == 1) then
        if (ltrans) then
            do i = 1, size(indices)
                j = indices(i)
                dst(:,j) = src(i,:)
            end do
        else
            do i = 1, size(indices)
                j = indices(i)
                dst(j,:) = src(i,:)
            end do
        end if
    else if (dim == 2) then
        if (ltrans) then
            do i = 1, size(indices)
                j = indices(i)
                dst(j,:) = src(:,i)
            end do
        else
            do i = 1, size(indices)
                j = indices(i)
                dst(:,j) = src(:,i)
            end do
        end if
    end if
end subroutine



subroutine unpack_indexed_1d (src, indices, dst, fill)
    !*  UNPACK_INDEXED unpacks an array into a (weakly) larger array
    !   using a given mapping of indices.
    !   Any element in the destination array that is not present in the
    !   source array is assigned a constant value, if present.
    real (PREC), intent(in), dimension(:), contiguous :: src
    integer, intent(in), dimension(:), contiguous :: indices
        !*  Mapping between SRC and DST. The value of SRC(i) will be
        !   copied to DST(INDICES(i)).
    real (PREC), intent(out), dimension(:) :: dst
    real (PREC), intent(in), optional :: fill
        !*  Fill value assigned to all elements in DST which are not present
        !   in SRC. If not present, missing elements will remain untouched.

    integer :: i, j

    if (present(fill)) then
        dst = fill
    end if

    do i = 1, size(indices)
        j = indices(i)
        dst(j) = src(i)
    end do

end subroutine




