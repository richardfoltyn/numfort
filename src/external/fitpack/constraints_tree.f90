module constraints_tree_mod
    implicit none
    private

    integer, parameter :: TREE_STATUS_SUCCESS = 0
    integer, parameter :: TREE_STATUS_FULL_CAPACITY = 1

    public :: TREE_STATUS_SUCCESS, TREE_STATUS_FULL_CAPACITY

    type, public :: constraints_tree
        private
        integer, dimension(:), pointer, contiguous :: up, left, right, info
        integer :: ifree = 1, ncount = 0, maxtr = -1, ibranch = -1
    contains
        procedure, public, pass :: init => tree_init
        procedure, public, pass :: remove_constraints
        procedure, public, pass :: add_constraint
        procedure, public, pass :: has_constraints
        procedure, public, pass :: next_constraint
        procedure, pass :: insert_node
        procedure, pass :: delete_node
        procedure, pass :: reset_branch_pointer
        procedure, pass :: remove_subtree
        procedure, pass :: insert_tail
    end type

contains

subroutine tree_init (self, maxtr, up, left, right, info)
    class (constraints_tree), intent(in out) :: self
    integer :: maxtr
    integer, dimension(maxtr), target :: up, left, right, info

    intent (in) :: maxtr, up, left, right, info

    integer :: status, idx

    self%up => up
    self%left => left
    self%right => right
    self%info => info
    self%maxtr = maxtr

    self%ifree = 1
    self%ncount = 0

    ! insert root node
    call insert_node (self, 1, 0, idx, status)

    ! set branch pointer to root node
    self%ibranch = 1

end subroutine

pure subroutine insert_node (self, parent, info, inew, status, leftof, rightof)
    class (constraints_tree), intent(in out) :: self
    integer, intent(in) :: parent, info, leftof, rightof
    integer, intent(out) :: status, inew

    optional :: leftof, rightof

    integer :: i, ifrom

    status = TREE_STATUS_SUCCESS

    if (self%ifree < 1 .or. self%ifree > self%maxtr) then
        status = TREE_STATUS_FULL_CAPACITY
        return
    end if

    inew = self%ifree

    self%up(inew) = parent
    self%left(inew) = 0
    self%right(inew) = 0
    self%info(inew) = info
    self%ncount = self%ncount + 1

    if (present(leftof)) then
        self%left(leftof) = inew
    else if (present(rightof)) then
        self%right(rightof) = inew
    end if

    ! find position of next free node
    ifrom = self%ifree + 1
    self%ifree = 0
    do i = ifrom, self%maxtr
        if (self%up(i) == 0) then
            self%ifree = i
            exit
        end if
    end do

end subroutine

pure subroutine delete_node (self, idx)
    class (constraints_tree), intent(in out) :: self
    integer, intent(in) :: idx

    ! do not allow root node to be deleted
    if (idx <= 1) return

    self%ncount = self%ncount - 1
    ! adjust index of first free node slot if necessary
    if (idx < self%ifree) self%ifree = idx

    self%up(idx) = 0
    self%left(idx) = 0
    self%right(idx) = 0
    self%info(idx) = 0

end subroutine

pure subroutine reset_branch_pointer (self)
    class (constraints_tree), intent(in out) :: self
    integer :: il

    ! point branch pointer to terminal node of left-most branch if there is one,
    ! and to root node otherwise
    self%ibranch = 1
    il = self%left(1)
    do while (il /= 0)
        self%ibranch = il
        il = self%left(il)
    end do
end subroutine

pure function has_constraints (self) result(res)
    class (constraints_tree), intent (in) :: self
    logical :: res

    res = (self%ibranch > 1)
end function

pure subroutine remove_constraints (self, nbind)
    class (constraints_tree), intent(in out) :: self
    integer, intent (in) :: nbind

    if (self%left(1) /= 0) then
        call remove_subtree(self, self%left(1), 1, nbind)
        ! remove reference to left child node if it was deleted
        if (self%up(self%left(1)) == 0) self%left(1) = 0
    end if

    call reset_branch_pointer (self)
end subroutine


pure recursive subroutine remove_subtree (self, idx, depth, nbind)
    class (constraints_tree), intent (in out) :: self
    integer, intent (in) :: idx, depth, nbind

    integer :: il, ir

    ir = self%right(idx)
    if (ir /= 0) then
        call remove_subtree (self, ir, depth, nbind)
        ! in case subtree on the right was removed, delete reference
        if (self%up(ir) == 0) self%right(idx) = 0
    end if

    il = self%left(idx)
    if (il /= 0) then
        call remove_subtree (self, il, depth + 1, nbind)
        if (self%up(il) == 0) self%left(idx) = 0
    end if

    ! if child (left) node was deleted, but right node is still valid
    ! and contains a sub-tree, then remove current node and put right node in
    ! its place.
    ! Don't do this at terminal nodes!
    il = self%left(idx)
    ir = self%right(idx)
    if (il == 0 .and. ir /= 0 .and. depth < nbind) then
        if (self%left(ir) /= 0 .or. self%right(ir) /= 0) then
            il = self%left(self%up(idx))
            ! left node of parent node, referenced directly
            if (il == idx) then
                self%left(self%up(idx)) = ir
            else
                ! need to locate left node which we have no direct
                ! reference to
                do while (self%right(il) /= idx)
                    il = self%right(il)
                end do
                ! remove current node from horizontal chain
                self%right(il) = self%right(idx)
            end if
            self%up(idx) = 0
            self%right(idx) = 0
        end if
    end if

    ! determine whether this is terminal node
    ! end of branch, no more nodes in subtree
    il = self%left(idx)
    ir = self%right(idx)
    if (il == 0 .and. ir == 0 .and. depth < nbind) then
        ! mark current node as deleted if required number of binding
        ! constraints not found
        call delete_node (self, idx)
    end if
end subroutine

pure subroutine add_constraint (self, ibind, status)
    class (constraints_tree), intent(in out) :: self
    integer, dimension(:), intent(in) :: ibind
    integer, intent(out) :: status

    integer :: il

    if (size(ibind) == 0) return

    il = self%left(1)
    if (il == 0) then
        call insert_node (self, 1, ibind(1), il, status, leftof=1)
        if (status /= TREE_STATUS_SUCCESS) return
    end if

    call insert_tail (self, il, ibind, status)

    call reset_branch_pointer (self)
end subroutine

recursive pure subroutine insert_tail (self, idx, ibind, status)
    class (constraints_tree), intent(in out) :: self
    integer :: idx, status
    integer, dimension(:) :: ibind

    intent (in) :: idx, ibind
    intent (out) :: status

    integer :: val, maxtr, il, inew

    val = self%info(idx)
    status = TREE_STATUS_SUCCESS

    if (val == ibind(1)) then
        if (size(ibind) > 1) then
            ! if child node exists, then pass processing of the remaining tail
            ! to child instance; otherwise create child node first
            if (self%left(idx) == 0) then
                ! need to create child note first
                call insert_node (self, idx, ibind(2), inew, status, leftof=idx)
                if (status /= TREE_STATUS_SUCCESS) goto 100
            end if

            call insert_tail (self, self%left(idx), ibind(2:), status)
            if (status /= TREE_STATUS_SUCCESS) goto 100
        ! else
        !     ! end of branch reached, nothing more to insert. Set current
        !     ! branch pointer to this index.
        !     self%ibranch = idx
        end if
    else if (val < ibind(1)) then
        ! move to the right, if necessary creating a right node first
        if (self%right(idx) == 0) then
            call insert_node (self, self%up(idx), ibind(1), inew, status, rightof=idx)
            if (status /= TREE_STATUS_SUCCESS) goto 100
        end if

        call insert_tail (self, self%right(idx), ibind, status)
        if (status /= TREE_STATUS_SUCCESS) goto 100
    else if (val > ibind(1)) then
        ! need to insert a new node before the current one
        ! first, find the node with right pointing to current node
        il = self%left(self%up(idx))
        if (il == idx) then
            ! current node is the first child
            call insert_node (self, self%up(idx), ibind(1), inew, status, leftof=self%up(idx))
        else
            ! Need to insert a sibling whose right reference points to current
            ! node. Note that even if there are other siblings to the right,
            ! it must be true that their info < val, as otherwise we wouldn't
            ! end up at this point.
            do while (self%right(il) /= idx)
                il = self%right(il)
            end do
            call insert_node (self, self%up(idx), ibind(1), inew, status, rightof=il)
        end if

        if (status /= TREE_STATUS_SUCCESS) goto 100
        self%right(inew) = idx

        call insert_tail (self, inew, ibind, status)
        if (status /= TREE_STATUS_SUCCESS) goto 100
    end if

    ! No special error handling required at this point
100 continue
end subroutine

pure subroutine next_constraint (self, ibind)
    class (constraints_tree), intent(in out) :: self
    integer, dimension(:), intent(in out) :: ibind

    integer :: inode, i, nbind, level
    
    ! branch length
    nbind = size(ibind)

    ! no valid pointer to last branch
    if (self%ibranch == 1 .or. nbind == 0) return
    
    inode = self%ibranch
    do i = nbind, 1, -1
        ibind(i) = self%info(inode)
        inode = self%up(inode)
    end do

    ! find new branch endpoint and store it in merk. We assume that *all*
    ! branches in tree are of length nbind, since fpseno is called right
    ! after fpdeno in fpcosp.
    inode = self%ibranch
    level = nbind
    do while (self%right(inode) == 0 .and. level > 0)
        inode = self%up(inode)
        level = level - 1
    end do

    if (inode == 1) then
        self%ibranch = 1
        return
    end if

    ! at this point we found a node with a sibling to the right
    inode = self%right(inode)
    do while (level < nbind)
        if (self%left(inode) /= 0) then
            inode = self%left(inode)
        else if (self%right(inode) /= 0) then
            inode  = self%right(inode)
        end if
        level = level + 1
    end do
    self%ibranch = inode

end subroutine

end module
