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
        procedure, public, pass :: keep_constraints
        procedure, public, pass :: add_constraint
        procedure, public, pass :: has_constraints
        procedure, public, pass :: next_constraint
        procedure, public, pass :: reset_branch_pointer
        procedure, pass :: insert_node
        procedure, pass :: delete_node
        procedure, pass :: remove_subtree
    end type

contains

subroutine tree_init (self, up, left, right, info)
    class (constraints_tree), intent(in out) :: self
    integer, dimension(:), target, contiguous :: up, left, right, info

    intent (in) :: up, left, right, info

    integer :: status, idx

    self%up => up
    self%left => left
    self%right => right
    self%info => info
    self%maxtr = size(up)

    self%ifree = 1
    self%ncount = 0

    ! erase remaining space to make sure there is no nodes from previous call
    ! of fpcosp
    self%up = 0
    self%left = 0
    self%right = 0
    self%info = 0

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

pure subroutine keep_constraints (self, nbind)
    class (constraints_tree), intent(in out) :: self
    integer, intent (in) :: nbind

    if (self%left(1) /= 0) then
        call remove_subtree(self, self%left(1), 1, nbind)
        ! remove reference to left child node if it was deleted
        if (self%up(self%left(1)) == 0) self%left(1) = 0
    end if
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

    integer :: nbind, inode, vnode, vinsert, inew, il, level

    status = TREE_STATUS_SUCCESS
    nbind = size(ibind)

    if (nbind == 0) return

    inode = self%left(1)
    level = 1

    if (inode == 0) then
        call insert_node (self, 1, ibind(1), inode, status, leftof=1)
        if (status /= TREE_STATUS_SUCCESS) return
        level = level + 1
    end if

    do while (level <= nbind)
        vnode = self%info(inode)
        vinsert = ibind(level)

        if (vnode == vinsert) then
            ! need to pass execution on to child node
            if (level < nbind) then
                ! if child node exists, then pass processing of the remaining tail
                ! to child instance; otherwise create child node first
                if (self%left(inode) == 0) then
                    ! need to create child note first
                    call insert_node (self, inode, vinsert, inew, status, leftof=inode)
                    if (status /= TREE_STATUS_SUCCESS) goto 100
                    inode = inew
                else
                    inode = self%left(inode)
                end if

                level = level + 1
            end if
        else if (vnode < vinsert) then
            ! move to the right, if necessary creating a right node first
            if (self%right(inode) == 0) then
                call insert_node (self, self%up(inode), vinsert, inew, status, rightof=inode)
                if (status /= TREE_STATUS_SUCCESS) goto 100
                inode = inew
            else
                inode = self%right(inode)
            end if
        else
            ! Handle case vnode > vinsert:
            ! need to insert a new node the current one
            ! first, so find the node with right pointing to current node
            il = self%left(self%up(inode))
            if (il == inode) then
                ! current node is the first child
                call insert_node (self, self%up(inode), vinsert, inew, status, leftof=self%up(inode))
            else
                ! Need to insert a sibling whose right reference points to current
                ! node. Note that even if there are other siblings to the right,
                ! it must be true that their info < val, as otherwise we wouldn't
                ! end up at this point.
                do while (self%right(il) /= inode)
                    il = self%right(il)
                end do
                call insert_node (self, self%up(inode), vinsert, inew, status, rightof=il)
            end if

            if (status /= TREE_STATUS_SUCCESS) goto 100

            self%right(inew) = inode
            inode = inew
        end if
    end do

    ! No special error handling required at this point
    return
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

    ! find new branch endpoint and store it in merk. This algorithm only
    ! searches for candidate branches to the right of the current branch.
    ! This should be sufficient since reset_branch_pointer() is called
    ! in fpcosp prior to first calling this routine for a particular number
    ! of binding constraints. Thus when next_constraint() is called for the
    ! first time for a particular value of nbind, the branch pointer should
    ! point to the left-most branch.
    inode = self%ibranch
    level = nbind
    start_search : do while (.true.)
        do while (self%right(inode) == 0 .and. level > 0)
            inode = self%up(inode)
            level = level - 1
        end do

        ! reached root node, no other branch found
        if (inode == 1) exit

        ! at this point we found a node with a sibling to the right;
        ! move left/right until we reached desired level; if branch is not
        ! long enough we need to go back up and move to the right.
        inode = self%right(inode)
        do while (level < nbind)
            if (self%left(inode) /= 0) then
                inode = self%left(inode)
            else if (self%right(inode) /= 0) then
                inode  = self%right(inode)
            else
                ! dead end: have not reached desired level yet, but cannot
                ! go either down or right; need to move back up until we
                ! find a sibling to the right
                cycle start_search
            end if
            level = level + 1
        end do
        ! check that the node found is at level nbind, ie does not have
        ! child nodes; in that case we are done.
        if (level == nbind .and. self%left(inode) == 0) exit
    end do start_search

    ! this should never happen!
    if (inode /= 1 .and. level /= nbind) then
        self%ibranch = 1
    end if

    self%ibranch = inode

end subroutine

end module
