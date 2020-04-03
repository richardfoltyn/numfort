

subroutine cast_to_args_default (tgt, ptr)
    !*  CAST_TO_ARGS_DEFAULT returns a pointer of type ARGS_DEFAULT if
    !   the given argument is of the respective type.
    class (args_data), intent(in), target :: tgt
    type (args_default), intent(inout), pointer :: ptr

    nullify (ptr)

    select type (tgt)
    type is (args_default)
        ptr => tgt
    end select

end subroutine


subroutine cond_alloc_args_default (self, nr, ni)
    !*  COND_ALLOC conditionally allocates the array attributes of
    !   and ARGS_DEFAULT object to the desired sizes.
    type (args_default), intent(inout) :: self
    integer, intent(in), optional :: nr
    integer, intent(in), optional :: ni

    if (present(nr)) then
        call cond_alloc (self%rdata, nr)
    end if

    if (present(ni)) then
        call cond_alloc (self%idata, nr)
    end if
end subroutine
