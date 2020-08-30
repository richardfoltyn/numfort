

module numfort_common_transpose

    use, intrinsic :: iso_fortran_env

    implicit none

    private

    public :: transpose_no_copy


    interface transpose_no_copy
        procedure transpose_no_copy_logical, &
            transpose_no_copy_real32, transpose_no_copy_real64
    end interface

    contains



subroutine transpose_no_copy_logical (src, dst)
    !*  TRANSPOSE_NO_COPY stores the transposed SRC to DST.
    !
    !   Routine is intended to avoid temporary arrays created by compiler
    !   if SRC has a TARGET attribute and DST is a pointer.
    !
    !   Assumes that SRC and DST do not point to same memory!
    logical, intent(in), dimension(:,:), contiguous :: src
    logical, intent(out), dimension(:,:), contiguous :: dst

    dst = transpose (src)
end subroutine



subroutine transpose_no_copy_real32 (src, dst)
    !*  TRANSPOSE_NO_COPY stores the transposed SRC to DST.
    !
    !   Routine is intended to avoid temporary arrays created by compiler
    !   if SRC has a TARGET attribute and DST is a pointer.
    !
    !   Assumes that SRC and DST do not point to same memory!
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:), contiguous :: src
    real (PREC), intent(out), dimension(:,:), contiguous :: dst

    dst = transpose (src)
end subroutine


subroutine transpose_no_copy_real64 (src, dst)
    !*  TRANSPOSE_NO_COPY stores the transposed SRC to DST.
    !
    !   Routine is intended to avoid temporary arrays created by compiler
    !   if SRC has a TARGET attribute and DST is a pointer.
    !
    !   Assumes that SRC and DST do not point to same memory!
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:), contiguous :: src
    real (PREC), intent(out), dimension(:,:), contiguous :: dst

    dst = transpose (src)
end subroutine

end module
