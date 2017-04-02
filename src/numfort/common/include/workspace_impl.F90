pure subroutine __APPEND(ws_assert_alloc,__PREC) (self, nrwrk, niwrk, ncwrk, nlwrk)

    integer, parameter :: PREC = __PREC
    ! Note: avoid polymorphic calls, might break OpenMP
    type (__APPEND(workspace,__PREC)), intent(in out) :: self
    integer, intent(in), optional :: nrwrk, niwrk, ncwrk, nlwrk

    real (PREC), dimension(:), allocatable :: rtmp
    integer, dimension(:), allocatable :: itmp
    logical, dimension(:), allocatable :: ltmp
    character (:), allocatable :: ctmp

    ! Real working array
    if (present(nrwrk)) then
        if (nrwrk > 0) then
            if (.not. allocated(self%rwrk)) then
                allocate (self%rwrk(nrwrk))
            else if (size(self%rwrk) < nrwrk) then
                allocate (rtmp(nrwrk))
                rtmp(1:size(self%rwrk)) = self%rwrk
                call move_alloc (rtmp, self%rwrk)
            end if
            self%nrwrk = size(self%rwrk)
        end if
    end if

    ! Integer working array
    if (present(niwrk)) then
        if (niwrk > 0) then
            if (.not. allocated(self%iwrk)) then
                allocate (self%iwrk(niwrk))
            else if (size(self%iwrk) < niwrk) then
                allocate (itmp(niwrk))
                itmp(1:size(self%iwrk)) = self%iwrk
                call move_alloc (itmp, self%iwrk)
            end if
            self%niwrk = size(self%iwrk)
        end if
    end if

    ! Logical working array
    if (present(nlwrk)) then
        if (nlwrk > 0) then
            if (.not. allocated(self%lwrk)) then
                allocate (self%lwrk(nlwrk))
            else if (size(self%lwrk) < nlwrk) then
                allocate (ltmp(nlwrk))
                ltmp(1:size(self%lwrk)) = self%lwrk
                call move_alloc (ltmp, self%lwrk)
            end if
            self%nlwrk = size(self%lwrk)
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
            self%ncwrk = len(self%cwrk)
        end if
    end if
end subroutine


pure subroutine __APPEND(ws_assert_alloc_ptr,__PREC) (ws, ptr_ws)
    !*  ASSERT_ALLOC_PTR ensures that ptr_ws points to allocated memory:
    !   If argument ws is present, ptr_ws points to ws and no additional
    !   memory is allocated. If ws is not present, then ptr_ws points to a newly
    !   allocated instance of ws.
    type (__APPEND(workspace,__PREC)), intent(in out), target, optional :: ws
    type (__APPEND(workspace,__PREC)), pointer, intent(in out) :: ptr_ws

    if (present(ws)) then
        ptr_ws => ws
    else
        allocate (ptr_ws)
    end if

end subroutine

pure subroutine __APPEND(ws_assert_dealloc_ptr,__PREC) (ws, ptr_ws)
    !*  ASSERT_DEALLOC_PTR ensures that memory that was dynamically allocated
    !   by ASSERT_ALLOC_PTR is released. The memory pointed to by ptr_ws
    !   needs to be released only if argument ws is not present, as only then
    !   a new array was allocated by ASSERT_ALLOC_PTR, which is referenced
    !   by ptr_ws.

    type (__APPEND(workspace,__PREC)), intent(in), target, optional :: ws
    type (__APPEND(workspace,__PREC)), pointer, intent(in out) :: ptr_ws

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