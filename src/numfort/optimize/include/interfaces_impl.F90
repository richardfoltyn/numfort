

subroutine __APPEND(cast_to_args_default,__PREC) (tgt, ptr)
    !*  CAST_TO_ARGS_DEFAULT returns a pointer of type ARGS_DEFAULT if
    !   the given argument is of the respective type.
    class (args_data), intent(in), target :: tgt
    type (__APPEND(args_default,__PREC)), intent(inout), pointer :: ptr

    nullify (ptr)

    select type (tgt)
    type is (__APPEND(args_default,__PREC))
        ptr => tgt
    end select

end subroutine


subroutine __APPEND(cond_alloc_args_default,__PREC) (self, nr, ni)
    !*  COND_ALLOC conditionally allocates the array attributes of
    !   and ARGS_DEFAULT object to the desired sizes.
    type (__APPEND(args_default,__PREC)), intent(inout) :: self
    integer, intent(in), optional :: nr
    integer, intent(in), optional :: ni

    if (present(nr)) then
        call cond_alloc (self%rdata, nr)
    end if

    if (present(ni)) then
        call cond_alloc (self%idata, nr)
    end if
end subroutine
