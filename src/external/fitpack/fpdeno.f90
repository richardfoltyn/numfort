subroutine fpdeno2 (maxtr, up, left, right, nbind, merk)
    implicit none
    integer :: nbind, merk, maxtr
    integer, dimension(maxtr) :: up, left, right

    intent (in) :: maxtr, nbind
    intent (in out) :: up, left, right, merk

    integer :: il

    if (left(1) /= 0) then
        if (up(left(1)) /= 0) then
            call check_subtree(left(1), 1, maxtr, up, left, right, nbind)
            ! remove reference to left child node if it was deleted
            ! if (up(left(1)) == 0) left(1) = 0
        end if
    end if

    ! point merk to terminal node of left-most branch if there is one, and to
    ! 1 otherwise
    merk = 1
    il = left(1)
    do while (il /= 0)
        merk = il
        il = left(il)
    end do

end subroutine

pure recursive subroutine check_subtree (idx, depth, maxtr, up, left, right, nbind)

    integer :: idx, depth, maxtr, nbind
    integer, dimension(maxtr) :: up, left, right

    intent (in) :: idx, depth, maxtr, nbind
    intent (in out) :: up, left, right

    integer :: il, ir
    logical :: delete, has_left, has_right

    ir = right(idx)
    if (ir /= 0) then
        ! make sure right node is valid (ie has up pointer)
        if (up(ir) /= 0) then
            call check_subtree (ir, depth, maxtr, up, left, right, nbind)
            ! in case subtree on the right was removed, delete reference
            ! if (up(ir) == 0) right(idx) = 0
        end if
    end if

    il = left(idx)
    if (il /= 0) then
        if (up(il) /= 0) then
            call check_subtree (il, depth + 1, maxtr, up, left, right, nbind)
            ! if (up(il) == 0) left(idx) = 0
        end if
    end if

    ! if child (left) node was deleted, but right node is still valid
    ! and contains a sub-tree, then
    ! remove current node and put right node in its place.
    ! Don't do this at terminal nodes!
    ! delete = .false.
    ! if (left(idx) == 0 .and. right(idx) /= 0 .and. depth < nbind) then
    !     ir = right(idx)
    !     if (up(ir) /= 0) then
    !         ! check whether right node has a valid left child
    !         if (left(ir) /= 0) then
    !             if (up(left(ir)) /= 0) delete = .true.
    !         end if
    !         ! check whether right node has a valid right child
    !         if (right(ir) /= 0) then
    !             if (up(right(ir)) /= 0) delete = .true.
    !         end if
    !         if (delete) then
    !             il = left(up(idx))
    !             ! left node of parent node, referenced directly
    !             if (il == idx) then
    !                 left(up(idx)) = ir
    !             else
    !                 ! need to locate left node which we have no direct
    !                 ! reference to
    !                 do while (right(il) /= idx)
    !                     il = right(il)
    !                 end do
    !                 ! remove current node from horizontal chain
    !                 right(il) = right(idx)
    !             end if
    !             up(idx) = 0
    !             right(idx) = 0
    !         end if
    !     end if
    ! end if

    ! determine whether this is terminal node
    has_left = (left(idx) /= 0)
    if (has_left) has_left = (up(left(idx)) /= 0)
    has_right = (right(idx) /= 0)
    if (has_right) has_right = (up(right(idx)) /= 0)

    ! end of branch, no more nodes in subtree
    if (.not. has_left .and.  .not. has_right) then
        ! mark current node as deleted if required number of binding
        ! constraints not found
        if (depth < nbind) up(idx) = 0
    end if
end subroutine
