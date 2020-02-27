


pure subroutine __APPEND(check_input_ext,__PREC) (rkind, ext, left, right, status)
    !*  CHECK_INPUT_EXT checks the `ext` argument that determines how
    !   values are extrapolated for many interpolation functions.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in) :: rkind
    integer (NF_ENUM_KIND), intent(in), optional :: ext
    real (PREC), intent(in), optional :: left, right
    type (status_t), intent(out) :: status

    status = NF_STATUS_OK

    if (present(ext)) then
        if (ext == NF_INTERP_EVAL_CONST) then
            if (.not. present(left) .or. .not. present(right)) then
                status = NF_STATUS_INVALID_ARG
                goto 100
            end if
        end if
    end if

100 continue

end subroutine
