
pure subroutine ws_reset (self)
    !*  WORKSPACE_CLEAR clear any internal state other than allocated
    !   working arrays.

    type (workspace), intent(inout) :: self

    self%roffset = 0
    self%ioffset = 0
    self%coffset = 0
    self%loffset = 0

end subroutine


pure subroutine ws_assert_alloc_int64 (self, nrwrk, niwrk, ncwrk, nlwrk)
    integer, parameter :: INTSIZE = int64
    type (workspace), intent(inout) :: self
    integer (INTSIZE), intent(in), optional :: nrwrk, niwrk, ncwrk, nlwrk

    real (PREC), dimension(:), allocatable :: rtmp
    integer, dimension(:), allocatable :: itmp
    logical, dimension(:), allocatable :: ltmp
    character (:), allocatable :: ctmp

    ! Real working array
    if (present(nrwrk)) then
        if (nrwrk > 0) then
            if (.not. allocated(self%rwrk)) then
                allocate (self%rwrk(nrwrk), source=0.0_PREC)
            else if (size(self%rwrk) < nrwrk) then
                allocate (rtmp(nrwrk), source=0.0_PREC)
                rtmp(1:size(self%rwrk)) = self%rwrk
                call move_alloc (rtmp, self%rwrk)
            end if
        end if
    end if

    ! Integer working array
    if (present(niwrk)) then
        if (niwrk > 0) then
            if (.not. allocated(self%iwrk)) then
                allocate (self%iwrk(niwrk), source=0)
            else if (size(self%iwrk) < niwrk) then
                allocate (itmp(niwrk), source=0)
                itmp(1:size(self%iwrk)) = self%iwrk
                call move_alloc (itmp, self%iwrk)
            end if
        end if
    end if

    ! Logical working array
    if (present(nlwrk)) then
        if (nlwrk > 0) then
            if (.not. allocated(self%lwrk)) then
                allocate (self%lwrk(nlwrk), source=.false.)
            else if (size(self%lwrk) < nlwrk) then
                allocate (ltmp(nlwrk), source=.false.)
                ltmp(1:size(self%lwrk)) = self%lwrk
                call move_alloc (ltmp, self%lwrk)
            end if
        end if
    end if

    if (present(ncwrk)) then
        if (ncwrk > 0) then
            if (.not. allocated(self%cwrk)) then
                allocate (character (ncwrk) :: self%cwrk)
            else if (len(self%cwrk) < ncwrk) then
                allocate (character (ncwrk) :: ctmp)
                call move_alloc (ctmp, self%cwrk)
            end if
        end if
    end if
end subroutine


pure subroutine ws_assert_alloc_int32 (self, nrwrk, niwrk, ncwrk, nlwrk)
    integer, parameter :: INTSIZE = int32
    type (workspace), intent(inout) :: self
    integer (INTSIZE), intent(in) :: nrwrk
    integer (INTSIZE), intent(in), optional :: niwrk, ncwrk, nlwrk

    integer (int64) :: lnrwrk, lniwrk, lncwrk, lnlwrk

    lnrwrk = int(nrwrk, int64)
    lniwrk = 0
    lncwrk = 0
    lnlwrk = 0

    if (present(niwrk)) lniwrk = int(niwrk, int64)
    if (present(ncwrk)) lncwrk = int(ncwrk, int64)
    if (present(nlwrk)) lnlwrk = int(nlwrk, int64)

    call assert_alloc (self, lnrwrk, lniwrk, lncwrk, lnlwrk)
end subroutine



pure subroutine ws_assert_alloc_ptr (ws, ptr_ws)
    !*  ASSERT_ALLOC_PTR ensures that ptr_ws points to allocated memory:
    !   If argument ws is present, ptr_ws points to ws and no additional
    !   memory is allocated. If ws is not present, then ptr_ws points to a newly
    !   allocated instance of ws.
    type (workspace), intent(inout), target, optional :: ws
    type (workspace), pointer, intent(inout) :: ptr_ws

    if (present(ws)) then
        ptr_ws => ws
    else
        allocate (ptr_ws)
    end if

end subroutine



pure subroutine ws_assert_dealloc_ptr (ws, ptr_ws)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_ws
    !   needs to be released only if argument ws is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_ws.
    type (workspace), intent(in), target, optional :: ws
    type (workspace), pointer, intent(inout) :: ptr_ws

    if (associated(ptr_ws)) then
        if (.not. present(ws)) then
            deallocate (ptr_ws)
            nullify (ptr_ws)
        else if (.not. associated(ptr_ws, ws)) then
            deallocate (ptr_ws)
            nullify (ptr_ws)
        end if
    end if
end subroutine



pure subroutine ws_get_rptr_1d_int32 (self, n, ptr, status)
    !*  WORKSPACE_GET_PTR associates a given pointer to a 1d-array with
    !   a subset of the workspace array of desired length.
    !
    !   Note: Cannot be implemented as function returning 1d-pointer,
    !   as it won't compile with ifort (CONTIGUOUS attribute on function
    !   return values seems to be ignored).

    integer, parameter :: INTSIZE = int32

    type (workspace), intent(inout), target :: self
    integer (INTSIZE), intent(in) :: n
    real (PREC), intent(out), dimension(:), pointer, contiguous :: ptr
    type (status_t), intent(out), optional :: status

    integer (int64) :: ifrom, ito
    type (status_t) :: lstatus

    nullify (ptr)
    lstatus = NF_STATUS_OK

    ifrom = self%roffset + 1
    ito = self%roffset + n

    ! Note: if working array is not large enough, return NULL pointer instead
    ! of reallocating the working array, which will invalidate all previosly
    ! associated pointers.
    if (ito > size(self%rwrk)) then
        lstatus = NF_STATUS_BOUNDS_ERROR
        goto 100
    end if

    ptr => self%rwrk(ifrom:ito)

    self%roffset = ito

100 continue
    if (present(status)) status = lstatus

end subroutine


pure subroutine ws_get_rptr_2d_int32 (self, shp, ptr, status)
    !*  WORKSPACE_GET_PTR maps a pointer to a 2d-array of given shape to
    !   a subset of the 1-dimensional workspace array of appropriate length.
    !
    !   Note: Cannot be implemented as function returning 1d-pointer,
    !   as it won't compile with ifort (CONTIGUOUS attribute on function
    !   return values seems to be ignored).
    integer, parameter :: INTSIZE = int32

    type (workspace), intent(inout), target :: self
    integer (INTSIZE), intent(in), dimension(:) :: shp
    real (PREC), intent(out), dimension(:,:), pointer, contiguous :: ptr
    type (status_t), intent(out), optional :: status

    integer (int64) :: ifrom, ito, n
    type (status_t) :: lstatus

    nullify (ptr)
    lstatus = NF_STATUS_OK

    n = product(shp)
    ifrom = self%roffset + 1
    ito = self%roffset + n

    ! Note: if working array is not large enough, return NULL pointer instead
    ! of reallocating the working array, which will invalidate all previosly
    ! associated pointers.
    if (ito > size(self%rwrk)) then
        lstatus = NF_STATUS_BOUNDS_ERROR
        goto 100
    end if

    ptr(1:shp(1),1:shp(2)) => self%rwrk(ifrom:ito)

    self%roffset = ito

100 continue
    if (present(status)) status = lstatus

end subroutine


pure subroutine ws_get_iptr_1d_int32 (self, n, ptr, status)
    !*  WORKSPACE_GET_PTR associates a given pointer to a 1d-array with
    !   a subset of the workspace array of desired length.
    !
    !   Note: Cannot be implemented as function returning 1d-pointer,
    !   as it won't compile with ifort (CONTIGUOUS attribute on function
    !   return values seems to be ignored).

    integer, parameter :: INTSIZE = int32

    type (workspace), intent(inout), target :: self
    integer (INTSIZE), intent(in) :: n
    integer, intent(out), dimension(:), pointer, contiguous :: ptr
    type (status_t), intent(out), optional :: status

    integer (int64) :: ifrom, ito
    type (status_t) :: lstatus

    nullify (ptr)
    lstatus = NF_STATUS_OK

    ifrom = self%ioffset + 1
    ito = self%ioffset + n

    ! Note: if working array is not large enough, return NULL pointer instead
    ! of reallocating the working array, which will invalidate all previosly
    ! associated pointers.
    if (ito > size(self%iwrk)) then
        lstatus = NF_STATUS_BOUNDS_ERROR
        goto 100
    end if

    ptr => self%iwrk(ifrom:ito)

    self%ioffset = ito

100 continue
    if (present(status)) status = lstatus

end subroutine
